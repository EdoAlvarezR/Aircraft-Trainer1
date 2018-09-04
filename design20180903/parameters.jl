##############################################################################
#   PARAMETERS
##############################################################################
# ------ CONSTANTS
g = 9.82            # (m/s^2) gravity
rhoinf= 1.225       # (kg/m^3) air density at sea level and 15Cdegs
mu = 1.846/10^5     # (kg/m*s) air dynamic viscosity

# ------ DESIGN PARAMETERS
# magVinf = 5.40    # (m/s) optimum cruise speed

# ------ MATERIALS
rho_c = 0.483       # (kg/m^2) Cloroplast area density
rho_f = 24.8        # (kg/m^3) HD EPS foam volumetric density
rho_spar = 50/1000/0.75     # (kg/m) Carbon fiber spar linear density

# ------ GEOMETRIC PARAMETERS
# # Fuselage
# w = 8/100         # (m) width
# l = (560.000-230.000)/1000          # (m) length

# Wing
b_w = 0.75          # (m) span
lambda_w = 0        # (deg) sweep
c_wroot = 0.20      # (m) root chord
t_wroot = 0.10      # (ratio) max thickness/chord at root
c_wtip = 0.10       # (m) tip chord
t_wtip = 0.20       # (ratio) max thickness/chord at root
# barc_w = (c_wtip+c_wroot)/2
# AR_w = b_w/c_wtip # Aspect ratio
# tr_w = c_wtip/c_wroot             # Taper ratio
twist_wroot = 0     # (deg) root twist
twist_wtip = -4     # (deg) tip twist
gamma_w = 5         # (deg) dihedral


# Horizontal stabilizer
b_h = b_w/4         # (m) span
lambda_h = 0        # (deg) sweep
c_hroot = b_h/5     # (m) root chord
t_hroot = 0.10      # (ratio) max thickness/chord at root
c_htip = c_hroot    # (m) tip chord
t_htip = t_hroot    # (ratio) max thickness/chord at root
# barc_h = (c_htip+c_hroot)/2
# AR_h = b_h/c_htip # Aspect ratio
# tr_h = c_htip/c_hroot             # Taper ratio
twist_hroot = 0     # (deg) root twist
twist_htip = 0      # (deg) tip twist
gamma_h = 0         # (deg) dihedral

# Vertical stabilizer
h_v = 0.75*b_h      # (m) height
lambda_v = 20       # (deg) sweep
c_vroot = 2*h_v/5   # (m) root chord
t_vroot = 0.10      # (ratio) max thickness/chord at root
c_vtip = c_vroot    # (m) tip chord
t_vtip = t_vroot    # (ratio) max thickness/chord at root
# barc_v = (c_vtip+c_vroot)/2
# AR_v = b_v/c_vtip # Aspect ratio
# tr_v = c_vtip/c_vroot             # Taper ratio
twist_vroot = 0     # (deg) root twist
twist_vtip = 0      # (deg) tip twist
gamma_v = 0         # (deg) dihedral


# # Center wing is same width as canard
# # centerwing = vlm.simpleWing(b_c, AR_w, tr_w, twist_wroot*180/pi,
# # lambda_w*180/pi, gamma_w*180/pi;
# # twist_tip=twist_wmid*180/pi, n=n_w, r=r_w)
# # Wing tips
# t_y_tip = b_w/2
# t_x_tip = b_w/2*tan(lambda_w) #t_y_tip*tan(lambda3)
# t_z_tip = 0.0#b_c/2*tan(lambda_w)+b_w/2*tan(lambda_w)
# t_c_tip = c_wtip
# t_twist_tip = twist_wtip
# t_y_mid = b_c/2 #starts at the center section
# t_x_mid = b_c/2*tan(lambda_w)
# t_z_mid = 0.0#b_c/2*tan(lambda_w*180/pi)
# t_c_mid = c_wmid
# t_twist_mid = twist_wmid
#
# Sref = 2*b_w*barc_w      # Reference area
#
# # ------ AIRPLANE ASSEMBLY
# x_w, y_w, z_w = l*2/4, 0, -w*1/8   # Wing's position in fuselage
# x_c, y_c, z_c = l*1/16, 0, z_w+vertical_offset  # Canard's position in fuselage
# # x_t, y_t, z_t = l-l_troot, 0, w/2  # Tail's position in fuselage
#
# # ------ USEFUL FUNCTIONS
# calc_Re(Vinf, l) = rhoinf*Vinf*l/mu
# calc_qinf(Vinf) = 1/2*rhoinf*Vinf.^2
# Vinfmin, Vinfmax = 3, 20 # (m/s)
#
#
# # ------ ASSUMPTIONS
# CLratio = 0.1    # Canard-wing distribution of lift
# #   CLratio = CLcanard/CLtot
#
# # Mass of electronic components
# M_base = (9 + 16 + 47 + 14 + 36 + 3 + 32 + 15 + 9*4 + 50)/1000
# # Mass of structural components
# M_f = rho_c * 4*w*l
# M_w = rho_f * b_w/2/cos(lambda_w) * barc_w^2 * t_w
# M_c = rho_f * b_c/2/cos(lambda_c) * barc_c^2 * t_c
# M_wl = rho_c * barc_wl * h_wl
# # M_t = rho_f * barc_t * h_t*t_t
# # Total mass
# Mtot = sum([M_base, M_f, M_w, M_c, M_wl])
# # Required lift in cruise
# L = Mtot*g
#
# println("Total mass: $(Mtot) (kg)")

####### END OF PARAMETERS ####################################################
