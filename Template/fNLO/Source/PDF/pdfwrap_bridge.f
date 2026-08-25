c     The Python exporter writes the LHAPDF 6 implementation with this external
c     ABI.  The bundled NNPDF backend uses this adapter to expose the same
c     selectable interface without putting generated includes in a module.
      subroutine pdfwrap
      use pdfwrap_module, only: configure_pdf
      implicit none
      include 'pdf.inc'
      include '../alfas.inc'

      call configure_pdf(pdlabel,asmz,nloop)
      end


      subroutine numberPDFm(idummy)
      use pdfwrap_module, only: numberpdfm_impl
      implicit none
      integer idummy

      call numberpdfm_impl(idummy)
      end


      subroutine initPDFm(idummy1,idummy2)
      use pdfwrap_module, only: initpdfm_impl
      implicit none
      integer idummy1,idummy2

      call initpdfm_impl(idummy1,idummy2)
      end


      subroutine initPDFsetbynamem(idummy,cdummy)
      use pdfwrap_module, only: initpdfsetbynamem_impl
      implicit none
      integer idummy
      character(len=*) cdummy

      call initpdfsetbynamem_impl(idummy,cdummy)
      end
