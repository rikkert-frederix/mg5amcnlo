c Legacy analysis ABI for the module-owned dummy implementation.

      subroutine analysis_begin(nwgt,weights_info)
      use analysis_dummy_module, only:
     &     module_analysis_begin => analysis_begin
      implicit none
      integer nwgt
      character*(*) weights_info(*)

      call module_analysis_begin()
      return
      end


      subroutine analysis_end()
      use analysis_dummy_module, only:
     &     module_analysis_end => analysis_end
      implicit none
      call module_analysis_end()
      return
      end


      subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
      use process_dimensions, only: event_capacity
      use analysis_dummy_module, only:
     &     module_analysis_fill => analysis_fill
      implicit none
      double precision p(0:4,event_capacity),wgts(*)
      integer istatus(event_capacity),ipdg(event_capacity),ibody

      call module_analysis_fill()
      return
      end
