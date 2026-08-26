      INTEGER IPOS, JPOS
      INTEGER FKS_I_D(%(nconfs)d), FKS_J_D(%(nconfs)d)
      INTEGER FKS_J_FROM_I_D(%(nconfs)d, NEXTERNAL, 0:NEXTERNAL)
      INTEGER PARTICLE_TYPE_D(%(nconfs)d, NEXTERNAL)
      INTEGER PDG_TYPE_D(%(nconfs)d, NEXTERNAL)
      LOGICAL NEED_COLOR_LINKS_D(%(nconfs)d)

data fks_i_D / %(fks_i_values)s /
data fks_j_D / %(fks_j_values)s /

data need_color_links_d / %(need_color_links)s /

%(fks_j_from_i_lines)s

C     Colour representation: octet = 8, triplet = 3, singlet = 1
%(col_lines)s

C     Particle type according to PDG
%(pdg_lines)s
