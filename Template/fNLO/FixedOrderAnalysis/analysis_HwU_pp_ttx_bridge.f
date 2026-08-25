c External analysis ABI for analysis_HwU_pp_ttx.f90.  Keep this fixed-form
c bridge with the selected module; only one analysis bridge may be linked.

      subroutine analysis_begin(nwgt,weights_info)
      use analysis_hwu_pp_ttx_module,
     &     only: module_analysis_begin => analysis_begin
      implicit none
      integer nwgt
      character(len=*) weights_info(*)
      call module_analysis_begin(nwgt,weights_info)
      end


      subroutine analysis_end(dummy)
      use analysis_hwu_pp_ttx_module,
     &     only: module_analysis_end => analysis_end
      implicit none
      double precision dummy
      call module_analysis_end(dummy)
      end


      subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
      use process_dimensions, only: nexternal
      use analysis_hwu_pp_ttx_module,
     &     only: module_analysis_fill => analysis_fill
      implicit none
      double precision p(0:4,nexternal),wgts(*)
      integer istatus(nexternal),ipdg(nexternal),ibody
      call module_analysis_fill(p,istatus,ipdg,wgts,ibody)
      end
