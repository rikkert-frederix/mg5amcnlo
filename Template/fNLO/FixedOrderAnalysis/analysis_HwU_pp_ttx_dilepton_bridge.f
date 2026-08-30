c External analysis ABI for the comprehensive dileptonic ttbar analysis.

      subroutine analysis_begin(nwgt,weights_info)
      use analysis_hwu_pp_ttx_dilepton_module,
     &     only: module_analysis_begin => analysis_begin
      implicit none
      integer nwgt
      character(len=*) weights_info(*)
      call module_analysis_begin(nwgt,weights_info)
      end


      subroutine analysis_end()
      use analysis_hwu_pp_ttx_dilepton_module,
     &     only: module_analysis_end => analysis_end
      implicit none
      call module_analysis_end()
      end


      subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
      use process_dimensions, only: nexternal
      use analysis_hwu_pp_ttx_dilepton_module,
     &     only: module_analysis_fill => analysis_fill
      implicit none
      double precision p(0:4,nexternal),wgts(*)
      integer istatus(nexternal),ipdg(nexternal),ibody
      call module_analysis_fill(p,istatus,ipdg,wgts,ibody)
      end


      subroutine analysis_fill_multiplicative(p,nparticles,istatus,
     &     ipdg,wgts,ibody)
      use analysis_hwu_pp_ttx_dilepton_module,
     &     only: module_analysis_fill => analysis_fill
      implicit none
      integer nparticles,istatus(nparticles),ipdg(nparticles),ibody
      double precision p(0:4,nparticles),wgts(*)
      call module_analysis_fill(p,istatus,ipdg,wgts,ibody)
      end
