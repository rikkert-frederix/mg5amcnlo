c Legacy analysis ABI for the module-owned dummy implementation.

      subroutine analysis_begin(nwgt,weights_info)
      use analysis_dummy_module, only:
     &     module_analysis_begin => analysis_begin
      implicit none
      integer nwgt
      character*(*) weights_info(*)

      call module_analysis_begin(nwgt,weights_info)
      return
      end


      subroutine analysis_end(xnorm)
      use analysis_dummy_module, only:
     &     module_analysis_end => analysis_end
      implicit none
      double precision xnorm

      call module_analysis_end(xnorm)
      return
      end


      subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
      use process_dimensions, only: nexternal
      use analysis_dummy_module, only:
     &     module_analysis_fill => analysis_fill
      implicit none
      double precision p(0:4,nexternal),wgts(*)
      integer istatus(nexternal),ipdg(nexternal),ibody

      call module_analysis_fill(p,istatus,ipdg,wgts,ibody)
      return
      end
