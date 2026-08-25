      subroutine init_fks_metadata_bridge()
      use fks_metadata, only: initialize_fks_metadata
      implicit none
      include 'nexternal.inc'
      include 'fks_info.inc'

c Copy only immutable generated lookup data.  Mutable FKS COMMON blocks
c remain owned by their existing chooser and computation routines.
      call initialize_fks_metadata(fks_i_d,fks_j_d,extra_cnt_d,
     &     isplitorder_born_d,isplitorder_cnt_d,fks_j_from_i_d,
     &     particle_type_d,pdg_type_d,split_type_d,
     &     need_color_links_d)
      return
      end
