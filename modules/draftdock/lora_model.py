import torch
from torch import nn
from boltz.model.modules.diffusionv2 import DiffusionModule


class LoRALayer(nn.Module):
    """LoRA layer for low-rank adaptation."""

    def __init__(self, base_layer: nn.Linear, rank: int):
        super().__init__()
        self.base_layer = base_layer
        self.lora_A = nn.Parameter(torch.randn(base_layer.in_features, rank) * 0.01)
        self.lora_B = nn.Parameter(torch.randn(rank, base_layer.out_features) * 0.01)

    def forward(self, x):
        # Apply the base layer and the LoRA adaptation
        base_output = self.base_layer(x)
        lora_output = x @ self.lora_A @ self.lora_B
        return base_output + lora_output


class LoRADiffusion(nn.Module):
    """LoRA fine-tuning for the Boltz diffusion model.
    should inheret the DiffusionModule class from boltz.model.modules.diffusionv2
    with related funcitons for the foward pass etc."""

    def __init__(self, base_model: DiffusionModule, lora_rank: int = 4):
        super().__init__()
        self.base_model = base_model
        self.lora_layers = nn.ModuleDict()

        # Inject LoRA layers into the linear layers of the base model
        for name, module in base_model.named_modules():
            if isinstance(module, nn.Linear):
                self.lora_layers[name] = LoRALayer(module, lora_rank)

    def forward(self, *args, **kwargs):
        # Forward pass through the base model with LoRA layers applied
        outputs = {}
        for name, module in self.base_model.named_modules():
            if name in self.lora_layers:
                outputs[name] = self.lora_layers[name](module, *args, **kwargs)
            else:
                outputs[name] = module(*args, **kwargs)
        return outputs

