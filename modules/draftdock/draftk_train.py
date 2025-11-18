import torch
from torch import nn, optim
from torch.utils.data import DataLoader
from lora_model import LoRADiffusion

boltz_path = os.getenv("BOLTZ_PATH", "/home/limcaoco/opt/boltz/src")
sys.path.append(boltz_path)

from boltz.model.modules.diffusionv2 import AtomDiffusion
from scoring_funcitons import center_of_mass_loss

def draftk_train(
    base_model: AtomDiffusion,
    train_dataloader: DataLoader,
    num_epochs: int,
    lora_rank: int = 4,
    T: int = 50,
    K: int = 10,
    learning_rate: float = 2e-4,
    weight_decay: float = 0.1,
    device: str = "cuda",
    wandb_project: str = "draftk_training",
    loss_function: str = "center_of_mass_loss",
    loss_function_kwargs: dict = None,
):
    # Initialize wandb
    import wandb
    wandb.init(project=wandb_project, config={
        "num_epochs": num_epochs,
        "lora_rank": lora_rank,
        "learning_rate": learning_rate,
        "weight_decay": weight_decay,
        "T": T,
        "K": K,
        "loss_function": loss_function,
    })

    # Initialize LoRA model
    lora_model = LoRADiffusion(base_model, lora_rank).to(device)
    optimizer = optim.AdamW(lora_model.parameters(), lr=learning_rate, weight_decay=weight_decay)

    # Freeze base model parameters
    for param in base_model.parameters():
        param.requires_grad = False

    # Define available loss functions
    loss_functions = {
        "center_of_mass_loss": center_of_mass_loss,
        # Add other loss functions here as needed
    }

    if loss_function not in loss_functions:
        raise ValueError(f"Loss function '{loss_function}' is not supported.")

    selected_loss_function = loss_functions[loss_function]
    loss_function_kwargs = loss_function_kwargs or {}

    # Training loop
    for epoch in range(num_epochs):
        lora_model.train()
        epoch_loss = 0.0

        for batch in train_dataloader:
            feats = batch["feats"].to(device)
            atom_mask = batch["atom_mask"].to(device)

            # Sample initial noise
            x_t = torch.randn_like(feats["coords"]).to(device)

            # Perform diffusion sampling over T steps
            sigmas = base_model.sample_schedule(T)
            for t in range(T):
                sigma_t = sigmas[t]
                with torch.no_grad() if t < T - K else torch.enable_grad():
                    x_t = x_t.detach() if t == T - K else x_t
                    x_t = lora_model.preconditioned_network_forward(
                        noised_atom_coords=x_t,
                        sigma=sigma_t,
                        network_condition_kwargs={
                            "feats": feats,
                            "atom_mask": atom_mask,
                        },
                    )

            # Decode final denoised coordinates
            final_denoised_coords = x_t

            # Compute reward using the selected loss function
            loss = selected_loss_function(final_denoised_coords, atom_mask, **loss_function_kwargs)

            # Backward pass and optimization
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            epoch_loss += loss.item()

        avg_loss = epoch_loss / len(train_dataloader)
        print(f"Epoch {epoch + 1}/{num_epochs}, Loss: {avg_loss}")

        # Log metrics to wandb
        wandb.log({"epoch": epoch + 1, "loss": avg_loss})

    print("Training complete.")
    wandb.finish()
