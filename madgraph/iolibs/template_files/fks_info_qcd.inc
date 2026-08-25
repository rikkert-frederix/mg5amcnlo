      INTEGER IPOS, JPOS
      INTEGER FKS_I_D(%(nconfs)d), FKS_J_D(%(nconfs)d)
      INTEGER EXTRA_CNT_D(%(nconfs)d)
      INTEGER ISPLITORDER_BORN_D(%(nconfs)d)
      INTEGER ISPLITORDER_CNT_D(%(nconfs)d)
      INTEGER FKS_J_FROM_I_D(%(nconfs)d, NEXTERNAL, 0:NEXTERNAL)
      INTEGER PARTICLE_TYPE_D(%(nconfs)d, NEXTERNAL)
      INTEGER PDG_TYPE_D(%(nconfs)d, NEXTERNAL)
      LOGICAL SPLIT_TYPE_D(%(nconfs)d, %(nsplitorders)d)
      LOGICAL NEED_COLOR_LINKS_D(%(nconfs)d)

data fks_i_D / %(fks_i_values)s /
data fks_j_D / %(fks_j_values)s /

data need_color_links_d / %(need_color_links)s /

%(fks_j_from_i_lines)s

C     Split type: the order of elements is %(splitorders_name)s
%(split_type_lines)s

C     Extra counterterm metadata. A zero identifier means that no extra
C     counterterm is needed.
data extra_cnt_D / %(extra_cnt_values)s /
data isplitorder_born_D / %(isplitorder_born_values)s /
data isplitorder_cnt_D / %(isplitorder_cnt_values)s /

C     Colour representation: octet = 8, triplet = 3, singlet = 1
%(col_lines)s

C     Particle type according to PDG
%(pdg_lines)s
