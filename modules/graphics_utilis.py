
def make_banner():
    T = '''
████████╗ 
╚══██╔══╝ 
   ██║    
   ██║    
   ██║    
   ╚═╝    
    '''
    R = '''
██████╗  
██╔══██╗ 
██████╔╝ 
██╔══██╗ 
██║  ██║ 
╚═╝  ╚═╝ 
    '''
    I = '''
████████╗ 
╚══██╔══╝ 
   ██║    
   ██║    
████████╗ 
╚═══════╝ 
    '''
    A = '''
 █████╗  
██╔══██╗ 
██╔══██║ 
██║  ██║ 
██║  ██║ 
╚═╝  ╚═╝ 
    '''
    G = '''
 ██████╗  
██╔════╝  
██║  ███╗ 
██║   ██║ 
╚██████╔╝ 
 ╚═════╝  
    '''
    E = '''
███████╗ 
██╔════╝ 
███████╗ 
██╔════╝ 
███████╗ 
╚══════╝ 
    '''
    D = '''
██████═╗ 
██║  ██║ 
██║  ██║ 
██║  ██║ 
██████╔╝ 
╚═════╝  
    '''

    T_lines = T.splitlines()
    R_lines = R.splitlines()
    I_lines = I.splitlines()
    A_lines = A.splitlines()
    G_lines = G.splitlines()
    E_lines = E.splitlines()
    D_lines = D.splitlines()

    # Combine the lines horizontally to spell "TRIAGED" with reduced spacing
    banner = "\n".join(
        t.ljust(8) + r.ljust(8) + i.ljust(8) + a.ljust(8) + g.ljust(8) + e.ljust(8) + d
        for t, r, i, a, g, e, d in zip(T_lines, R_lines, I_lines, A_lines, G_lines, E_lines, D_lines)
    )

    banner_foot = """
\033[1;35m█████████████████████████████████████████████████████████████████
\033[1;35m        TRIAGED: The All-in-One Virtual Screening Pipeline       
\033[1;35m  Targeted Ranking of In silico and Artificially GEnerated Drugs
\033[0;37m                          Version 1.0.0
\033[1;37m               MAOM Lab 2025 - University of Michigan
\033[1;35m█████████████████████████████████████████████████████████████████"""

    return banner, banner_foot