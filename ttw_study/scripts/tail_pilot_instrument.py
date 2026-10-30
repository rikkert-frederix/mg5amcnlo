"""Diagnostic hooks applied only to the isolated pilot export."""
from pathlib import Path
from campaign import digest

DECLARATIONS = '''
  ! Bounded diagnostic state: no random draws and no integrand changes.
  integer(kind=8), save :: tail_n=0_8, tail_positive=0_8, tail_negative=0_8
  integer, save :: tail_records=0, tail_bad=0, tail_unit=-1, tail_prod_ifks=0
  integer, save :: tail_leaf_count=0
  double precision, save :: tail_max=0d0, tail_sum=0d0, tail_sq=0d0, tail_abs=0d0
  double precision, save :: tail_prod(6)=0d0, tail_leaves(512)=0d0
  double precision, save :: tail_momenta(0:3,6)=0d0, tail_model(11)=0d0
  integer(kind=8), save :: tail_masks(512)=0_8
  double precision, save :: tail_component_sum(32)=0d0, tail_component_sq(32)=0d0
  double precision, save :: tail_component_abs(32)=0d0
  public :: tail_begin, tail_record, tail_flush
'''

PROCEDURES = '''
  subroutine tail_begin()
    tail_prod_ifks=0
    tail_prod=0d0
    tail_momenta=0d0
    tail_model=0d0
    tail_leaf_count=0
  end subroutine

  subroutine tail_production(ifks,volume)
    use factorized_phase_space, only: factorized_radiation_state, fetch_factorized_radiation_state, &
                                     fetch_factorized_block_momenta
    integer, intent(in) :: ifks
    double precision, intent(in) :: volume
    type(factorized_radiation_state) :: state
    logical :: available
    tail_prod_ifks=ifks
    tail_prod=0d0
    tail_prod(1)=1d0/volume
    call fetch_factorized_radiation_state(real_event,0,state,available)
    if (available) tail_prod(2:6)=[state%xi,1d0-state%y,state%xi_hat,state%xi_max,state%shat]
    call fetch_factorized_block_momenta(real_event,0,6,tail_momenta,available)
  end subroutine

  subroutine tail_record(x,volume,values)
    use weight_lines, only: itype, bundle_component, nFKS, wgts, wgt, bjx, scales2
    use ieee_arithmetic, only: ieee_is_finite
    double precision, intent(in) :: x(ndimmax),volume,values(nintegrals)
    double precision :: value, absolute, component
    integer :: i,c,j
    if (tail_unit==-1) open(newunit=tail_unit,file='record_breakers.dat',status='replace')
    tail_n=tail_n+1_8
    value=values(2)
    if (.not.ieee_is_finite(value)) then
      tail_bad=tail_bad+1
      if (tail_bad<=64) write(tail_unit,*) 'NONFINITE',tail_n,value
      return
    end if
    tail_sum=tail_sum+value
    tail_sq=tail_sq+value*value
    tail_abs=tail_abs+abs(value)
    if (value>0d0) tail_positive=tail_positive+1_8
    if (value<0d0) tail_negative=tail_negative+1_8
    do c=1,min(32,nintegrals-first_bundle_component_integral+1)
      component=values(first_bundle_component_integral+c-1)
      tail_component_sum(c)=tail_component_sum(c)+component
      tail_component_sq(c)=tail_component_sq(c)+component*component
      tail_component_abs(c)=tail_component_abs(c)+abs(component)
    end do
    if (abs(value)<=tail_max) return
    tail_max=abs(value)
    if (tail_records>=64) return
    tail_records=tail_records+1
    if (multiplicative_nlo_enabled()) then
      absolute=sum(abs(tail_leaves(1:tail_leaf_count)))
    else
      absolute=sum(abs(wgts(1,1:icontr)))
    end if
    write(tail_unit,'(a,4(1x,i0),2(1x,es25.17e3))') 'POINT',tail_n,ichan,iconfig,ndim,value,absolute
    write(tail_unit,'(a,1x,es25.17e3)') 'VOLUME',volume
    write(tail_unit,'(a,*(1x,es25.17e3))') 'X',x(1:ndim)
    write(tail_unit,'(a,1x,i0,6(1x,es25.17e3))') 'PRODUCTION',tail_prod_ifks,tail_prod
    write(tail_unit,'(a,11(1x,es25.17e3))') 'MODEL',tail_model
    do i=1,6
      write(tail_unit,'(a,1x,i0,4(1x,es25.17e3))') 'MOMENTA',i,tail_momenta(:,i)
    end do
    if (multiplicative_nlo_enabled()) then
      do i=1,tail_leaf_count
        write(tail_unit,'(a,1x,i0,1x,es25.17e3)') 'LEAF',tail_masks(i),tail_leaves(i)
      end do
    else
      do i=1,icontr
        write(tail_unit,'(a,4(1x,i0),8(1x,es25.17e3))') 'LINE',i,itype(i),bundle_component(i),nFKS(i), &
          wgts(1,i),wgt(1:3,i),bjx(1:2,i),scales2(2:3,i)
      end do
    end if
    write(tail_unit,'(a)') 'END'
    flush(tail_unit)
  end subroutine

  subroutine tail_flush()
    integer :: c
    if (tail_unit==-1) return
    write(tail_unit,'(a,1x,i0,3(1x,es25.17e3),3(1x,i0))') &
      'SUMMARY',tail_n,tail_sum,tail_sq,tail_abs,tail_positive,tail_negative,tail_bad
    do c=1,min(32,nintegrals-first_bundle_component_integral+1)
      if (tail_component_abs(c)==0d0) cycle
      write(tail_unit,'(a,2(1x,i0),3(1x,es25.17e3))') 'MOMENT',c,tail_n, &
        tail_component_sum(c),tail_component_sq(c),tail_component_abs(c)
    end do
    close(tail_unit)
  end subroutine
'''


def replace(text, old, new):
    if text.count(old)!=1:
        raise ValueError('Diagnostic hook is not unique: '+old[:80])
    return text.replace(old,new)


def instrument(process):
    process=Path(process)
    source=process/'SubProcesses/driver_mintFO.f90'
    text=source.read_text()
    text=replace(text,'  public :: run_mintfo_driver\n',DECLARATIONS+'  public :: run_mintfo_driver\n')
    text=replace(text,'  subroutine test_multiplicative_leaf_cuts( &',PROCEDURES+'\n  subroutine test_multiplicative_leaf_cuts( &')
    text=replace(text,'      call mint(sigint)\n','      call mint(sigint)\n      call tail_flush()\n')
    # Capture production before the following decay maps overwrite globals.
    old='''    if (p_born(0, 1) < 0d0) return

    call compute_prefactors_n1body(vegas_wgt)'''
    new='''    if (p_born(0, 1) < 0d0) return
    if (.not.contribution_is_nlo_decay(contribution_for_fks(ifks))) &
      call tail_production(ifks,volume)

    call compute_prefactors_n1body(vegas_wgt)'''
    text=replace(text,old,new)
    # Model values at the actual production-real coupling, before decay maps.
    start=text.index('  subroutine evaluate_additive_n1body_sector(')
    end=text.index('  end subroutine evaluate_additive_n1body_sector',start)
    part=text[start:end]
    part=replace(part,'      call compute_real_emission()',
        '      if (.not.contribution_is_nlo_decay(contribution_for_fks(ifks))) &\n' +
        '        call tail_model_values(tail_model)\n      call compute_real_emission()')
    text=text[:start]+part+text[end:]
    # Pi's real map is available before the local density is captured.
    old='''        if (.not. (production_contribution .and. &
                   reuse_folded_production)) then
          call compute_prefactors_n1body(vegas_wgt)'''
    new='''        if (production_contribution) &
          call tail_production(configurations(contribution),volumes(contribution))
        if (.not. (production_contribution .and. &
                   reuse_folded_production)) then
          call compute_prefactors_n1body(vegas_wgt)'''
    text=replace(text,old,new)
    start=text.index('  double precision function sigint_multiplicative_impl(')
    end=text.index('  end function sigint_multiplicative_impl',start)
    part=text[start:end]
    part=replace(part,'          call compute_real_emission()',
        '          if (production_contribution) call tail_model_values(tail_model)\n' +
        '          call compute_real_emission()')
    text=text[:start]+part+text[end:]
    text=replace(text,'      total_weights = total_weights + leaf_weights\n','''      total_weights = total_weights + leaf_weights
      if (tail_leaf_count>=512) error stop 'Diagnostic leaf capacity exceeded'
      tail_leaf_count=tail_leaf_count+1
      tail_leaves(tail_leaf_count)=leaf_weights(1)
      tail_masks(tail_leaf_count)=leaf_mask
''')
    source.write_text(text)
    bridge=process/'SubProcesses/driver_mintFO_bridge.f'
    text=bridge.read_text()
    text=replace(text,'      use driver_mintfo_module, only: sigint_impl\n',
        '      use driver_mintfo_module, only: sigint_impl,tail_begin,\n     $     tail_record\n')
    text=replace(text,'      sigint=sigint_impl(', '      call tail_begin()\n      sigint=sigint_impl(')
    text=replace(text,'     $     abrv,wgt_me_born,wgt_me_real)\n      return',
        '     $     abrv,wgt_me_born,wgt_me_real)\n      call tail_record(xx,vegas_wgt,f)\n      return')
    text+='''
      subroutine tail_model_values(values)
      implicit none
      include 'coupl.inc'
      double precision values(11)
      values(1)=MDL_MT
      values(2)=MDL_MW
      values(3)=MDL_WT
      values(4)=MDL_WW
      values(5)=G
      values(6)=dble(GC_4)
      values(7)=dimag(GC_4)
      values(8)=dble(GC_5)
      values(9)=dimag(GC_5)
      values(10)=dble(GC_11)
      values(11)=dimag(GC_11)
      end
'''
    bridge.write_text(text)
    weights=process/'SubProcesses/fks_weights.f90'
    text=weights.read_text()
    text=replace(text,'! Check for NaN\'s and INF\'s. Simply skip the contribution', '''! Diagnostic only: count raw nonfinite lines without altering native handling.
    call tail_nonfinite_notice(wgt1,wgt2,wgt3)
! Check for NaN's and INF's. Simply skip the contribution''')
    text+='''
subroutine tail_nonfinite_notice(a,b,c)
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
  double precision, intent(in) :: a,b,c
  integer, save :: n=0,u=-1
  if (ieee_is_finite(a).and.ieee_is_finite(b).and.ieee_is_finite(c)) return
  n=n+1
  if (n>64) return
  if (u==-1) open(newunit=u,file='nonfinite.dat',status='replace')
  write(u,*) 'NONFINITE',n,a,b,c
  flush(u)
end subroutine
'''
    weights.write_text(text)
    return dict(max_records_per_worker=64,maximum_leaf_values=512,
        sampled_points='All point moments; bounded record-large net weights only.',
        replay='Same seed, input grids and fixed-count worker; coordinates alone are insufficient.',
        files={str(p.relative_to(process)):digest(p) for p in (source,bridge,weights)})
