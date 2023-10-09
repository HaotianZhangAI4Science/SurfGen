try:
    from .generate_prot_ply import compute_inp_surface
except:
    from utils.generate_prot_ply import compute_inp_surface

__all__ = ['compute_inp_surface']