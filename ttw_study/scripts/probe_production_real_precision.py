#!/usr/bin/env python3
"""Compare production real spin densities in double and quadruple arithmetic.

The standalone probe is built in a new directory, links existing HELAS
libraries read-only, and does not alter a process executable or campaign.
It promotes the generated production tree routines, not the sampled momenta.
Thus it tests amplitude arithmetic, Ward identities and azimuthal invariance,
but not the accuracy of phase-space generation or subtraction counterterms.

Input format: point count, followed for every point by ID and real-matrix
index (1,2,3), a line MT MW WT WW G, six real numbers specifying real and
imaginary GC4, GC5 and GC11, and six four-vector lines (E,px,py,pz).
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess


DEFAULT_SUBPROCESS = 'P0_dxu_ttxwp_t_wpb_wp_mupvm_tx_wmbx_wm_emvex_wp_mupvm'
FIELDS = ('trace_dp', 'trace_qp', 'relative_trace_error',
          'relative_density_frobenius_error', 'relative_density_max_error',
          'ward_density_ratio_dp', 'ward_density_ratio_qp',
          'rotated_trace_relative_error_dp', 'rotated_trace_relative_error_qp',
          'max_shell_residual_GeV2', 'max_momentum_residual_GeV')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def promote(text):
    """Promote a pure generated density/JAMP routine and its HELAS calls."""
    text = re.sub(r'\bSDM_PRODUCTION_REAL_', 'MP_SDM_PRODUCTION_REAL_', text)
    text = re.sub(r'\bREAL\*8\b', 'REAL*16', text)
    text = re.sub(r'\bCOMPLEX\*16\b', 'COMPLEX*32', text)
    text = re.sub(r'\bDOUBLE PRECISION\b', 'REAL*16', text)
    text = text.replace("INCLUDE 'coupl.inc'", "INCLUDE 'mp_coupl_same_name.inc'")
    text = re.sub(r'(?i)\bDBLE\(([^()]*)\)', r'REAL(\1,KIND=16)', text)
    text = text.replace('DCONJG(', 'CONJG(')
    text = re.sub(r'(?i)(\d(?:[\d.]*))[d]([+-]?\d+)', r'\1Q\2', text)
    text = re.sub(r'\bCALL (?!MP_)([A-Z][A-Z0-9_]*)\(', r'CALL MP_\1(', text)
    return text


def ward_hook(text, index, quadruple=False):
    """Replace the external gluon polarization by p/E on demand."""
    leg = {1: 6, 2: 1, 3: 2}[index]
    prefix = 'MP_' if quadruple else ''
    marker = "      INCLUDE 'mp_coupl_same_name.inc'" if quadruple else "      INCLUDE 'coupl.inc'"
    if text.count(marker) != 1:
        raise ValueError('Expected a single coupling include in JAMP routine')
    text = text.replace(marker, '      INTEGER PROBE_GAUGE_MODE\n'
                        '      COMMON/PROBE_GAUGE/PROBE_GAUGE_MODE\n'+marker)
    pattern = rf'(?m)^(      CALL {prefix}VXXXXX\(P\(0,{leg}\),ZERO,.*)$'
    replacement = (r'\1\n      IF (PROBE_GAUGE_MODE.EQ.1) THEN\n'
                   f'        W(5:8,{leg})=P(0:3,{leg})/P(0,{leg})\n'
                   '      ENDIF')
    text, count = re.subn(pattern, replacement, text)
    if count != 1:
        raise ValueError('Expected a single external gluon wavefunction')
    return text


def driver_source():
    calls = '\n'.join(f'case({i})\n  call sdm_production_real_{i}(p,0,r)\n'
                      f'  call mp_sdm_production_real_{i}(q,0,s)' for i in range(1, 4))
    return '''program real_precision_probe
  implicit none
  include 'coupl.inc'
  include 'mp_coupl.inc'
  integer :: n,k,id,index,j,probe_gauge_mode
  common/probe_gauge/probe_gauge_mode
  real(8) :: p(0:3,6),rotated(0:3,6),v(6),td,normd,wd,rd,angle
  real(16) :: q(0:3,6),qr(0:3,6),tq,normq,wq,rq,err,errmax,shell,cons,masses(6)
  complex(8) :: r(2,12,12),r0(2,12,12)
  complex(16) :: s(2,12,12),s0(2,12,12)
  read(*,*) n
  do k=1,n
    read(*,*) id,index
    read(*,*) mdl_mt,mdl_mw,mdl_wt,mdl_ww,g
    read(*,*) v
    gc_4=cmplx(v(1),v(2),kind=8)
    gc_5=cmplx(v(3),v(4),kind=8)
    gc_11=cmplx(v(5),v(6),kind=8)
    do j=1,6
      read(*,*) p(:,j)
    enddo
    mp__mdl_mt=real(mdl_mt,16)
    mp__mdl_mw=real(mdl_mw,16)
    mp__mdl_wt=real(mdl_wt,16)
    mp__mdl_ww=real(mdl_ww,16)
    mp__g=real(g,16)
    mp__gc_4=cmplx(gc_4,kind=16)
    mp__gc_5=cmplx(gc_5,kind=16)
    mp__gc_11=cmplx(gc_11,kind=16)
    q=real(p,16)
    probe_gauge_mode=0
    call evaluate(index,p,q,r,s)
    r0=r
    s0=s
    td=0d0
    tq=0q0
    do j=1,12
      td=td+real(r(1,j,j),8)
      tq=tq+real(s(1,j,j),16)
    enddo
    normd=sqrt(sum(abs(r(1,:,:))**2))
    normq=sqrt(sum(abs(s(1,:,:))**2))
    err=sqrt(sum(abs(cmplx(r(1,:,:),kind=16)-s(1,:,:))**2))/max(normq,tiny(1q0))
    errmax=maxval(abs(cmplx(r(1,:,:),kind=16)-s(1,:,:)))/max(maxval(abs(s(1,:,:))),tiny(1q0))
    probe_gauge_mode=1
    call evaluate(index,p,q,r,s)
    wd=sqrt(sum(abs(r(1,:,:))**2))/max(normd,tiny(1d0))
    wq=sqrt(sum(abs(s(1,:,:))**2))/max(normq,tiny(1q0))
    probe_gauge_mode=0
    angle=.713123456789d0
    rotated=p
    rotated(1,:)=cos(angle)*p(1,:)-sin(angle)*p(2,:)
    rotated(2,:)=sin(angle)*p(1,:)+cos(angle)*p(2,:)
    qr=q
    qr(1,:)=cos(real(angle,16))*q(1,:)-sin(real(angle,16))*q(2,:)
    qr(2,:)=sin(real(angle,16))*q(1,:)+cos(real(angle,16))*q(2,:)
    call evaluate(index,rotated,qr,r,s)
    rd=0d0
    rq=0q0
    do j=1,12
      rd=rd+real(r(1,j,j),8)
      rq=rq+real(s(1,j,j),16)
    enddo
    masses=[0q0,0q0,mp__mdl_mt,mp__mdl_mt,mp__mdl_mw,0q0]
    shell=maxval(abs(q(0,:)**2-sum(q(1:3,:)**2,dim=1)-masses**2))
    cons=maxval(abs(q(:,1)+q(:,2)-sum(q(:,3:6),dim=2)))
    write(*,'(a,2(1x,i0),11(1x,es26.17e3))') 'RESULT',id,index,td,real(tq,8), &
      real(abs(real(td,16)-tq)/max(abs(tq),tiny(1q0)),8),real(err,8),real(errmax,8),wd,real(wq,8), &
      abs(rd-td)/max(abs(td),tiny(1d0)),real(abs(rq-tq)/max(abs(tq),tiny(1q0)),8), &
      real(shell,8),real(cons,8)
  enddo
contains
subroutine evaluate(index,p,q,r,s)
  integer,intent(in) :: index
  real(8),intent(in) :: p(0:3,6)
  real(16),intent(in) :: q(0:3,6)
  complex(8),intent(out) :: r(2,12,12)
  complex(16),intent(out) :: s(2,12,12)
  select case(index)
'''+calls+'''
case default
  error stop 'Unsupported production real index'
end select
end subroutine
end program
'''


def fixed_declarations_to_free(path):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line[0] in 'cC*!':
            continue
        if len(line) > 5 and line[5] != ' ':
            rows[-1] += line[6:]
        else:
            rows.append(line[6:])
    return '\n'.join(rows)


def build(args):
    process = args.process.resolve()
    directory = args.output_dir.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    source = process/'SubProcesses'/args.subprocess
    generated = []
    originals = {}
    for index in range(1, 4):
        original = source/f'spin_density_production_real_{index}.f'
        text = original.read_text()
        # A probe must evaluate every helicity, independently of history.
        text = text.replace('IF (GOODHEL(H).OR.NTRY.LE.2) THEN', 'IF (.TRUE.) THEN')
        original_copy = directory/f'production_real_{index}_dp.f'
        quad_copy = directory/f'production_real_{index}_qp.f'
        original_copy.write_text(ward_hook(text, index))
        quad_copy.write_text(ward_hook(promote(text), index, quadruple=True))
        generated.extend((original_copy, quad_copy))
        originals[str(original)] = digest(original)
    driver = directory/'real_precision_probe.f90'
    driver_text = driver_source()
    for name in ('coupl.inc', 'mp_coupl.inc'):
        driver_text = driver_text.replace("  include '"+name+"'",
            fixed_declarations_to_free(process/'Source/MODEL'/name))
    driver.write_text(driver_text)
    generated.append(driver)
    binary = directory/'real_precision_probe'
    command = [args.compiler, '-O', '-ffixed-line-length-132', '-ffree-line-length-none',
               '-I'+str(process/'Source/MODEL'), *(str(p) for p in generated),
               str(process/'lib/libdhelas.a'), '-o', str(binary)]
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (directory/'build.log').write_text(result.stdout)
    manifest = dict(schema='ttw_production_real_precision_probe_v1', process=str(process),
        subprocess=args.subprocess, command=command, returncode=result.returncode,
        source_sha256=originals, generated_sha256={str(p): digest(p) for p in generated},
        library_sha256=digest(process/'lib/libdhelas.a'), fields=FIELDS,
        scope='Amplitude arithmetic on identical double-precision input momenta; '
              'whole production spin density, Ward test, and azimuthal trace invariance. '
              'No test of the phase-space map or complete subtracted integrand.')
    (directory/'build.json').write_text(json.dumps(manifest, indent=2)+'\n')
    if result.returncode:
        raise RuntimeError('Build failed: '+str(directory/'build.log'))
    print(binary)


def run(args):
    if args.output.exists():
        raise ValueError('Output already exists')
    result = subprocess.run([str(args.binary.resolve())], input=args.input.read_text(),
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
    rows = []
    for line in result.stdout.splitlines():
        tokens = line.split()
        if not tokens or tokens[0] != 'RESULT':
            continue
        if len(tokens) != 3+len(FIELDS):
            raise ValueError('Malformed probe output: '+line)
        values = [float(x.replace('D', 'E')) for x in tokens[3:]]
        if not all(math.isfinite(x) for x in values):
            raise ValueError('Nonfinite probe output')
        rows.append(dict(point=int(tokens[1]), matrix_index=int(tokens[2]), **dict(zip(FIELDS, values))))
    count = int(args.input.read_text().split()[0])
    report = dict(schema='ttw_production_real_precision_result_v1', binary=str(args.binary.resolve()),
        binary_sha256=digest(args.binary), input=str(args.input.resolve()),
        input_sha256=digest(args.input), returncode=result.returncode, expected_points=count,
        points=rows, stdout=result.stdout,
        interpretation='Quadruple arithmetic evaluates the same generated tree formulas on '
            'promoted double-precision momenta/couplings. Density differences constrain '
            'production real amplitude arithmetic, not complete subtraction or map accuracy. '
            'Ward density ratios are squared-amplitude diagnostics, not amplitude ratios.')
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    if result.returncode or len(rows) != count:
        raise RuntimeError('Probe failed or incomplete: '+str(args.output))
    print(args.output)


def self_test_input(args):
    if args.output.exists():
        raise ValueError('Output already exists')
    vectors = [(100., 0., 50.), (-30., 70., -10.), (10., -20., 40.), (-80., -50., -80.)]
    masses = (172.5, 172.5, 80.4, 0.)
    lo, hi = 0., 10.
    for _ in range(100):
        mid = (lo+hi)/2
        energy = sum(math.sqrt(m*m+mid*mid*sum(x*x for x in v)) for m, v in zip(masses, vectors))
        if energy < 1000:
            lo = mid
        else:
            hi = mid
    scale = (lo+hi)/2
    momenta = [(500., 0., 0., 500.), (500., 0., 0., -500.)]
    momenta += [(math.sqrt(m*m+scale*scale*sum(x*x for x in v)), *(scale*x for x in v))
                for m, v in zip(masses, vectors)]
    g = math.sqrt(4*math.pi*.118)
    rows = ['3']
    for index in range(1, 4):
        rows += [f'{index} {index}', f'172.5 80.4 0 0 {g:.17e}',
                 f'{-g:.17e} 0 0 {g:.17e} 0 .46']
        rows += [' '.join(f'{x:.17e}' for x in p) for p in momenta]
    args.output.write_text('\n'.join(rows)+'\n')
    print(args.output)


def extract(args):
    """Select the largest captured net weights, without rerunning integration."""
    if args.output.exists() or args.output.with_suffix('.json').exists():
        raise ValueError('Output already exists')
    candidates = []
    for path in args.traces:
        current = None
        for line in path.read_text().splitlines():
            fields = line.split()
            if not fields:
                continue
            tag = fields[0]
            if tag == 'POINT':
                current = dict(source=str(path.resolve()), ordinal=int(fields[1]),
                               net=float(fields[5]), momenta={})
            elif tag == 'MODEL' and current is not None:
                current['model'] = [float(x.replace('D', 'E')) for x in fields[1:]]
            elif tag == 'MOMENTA' and current is not None:
                current['momenta'][int(fields[1])] = [float(x.replace('D', 'E')) for x in fields[2:]]
            elif tag == 'PRODUCTION' and current is not None:
                current['fks'] = int(fields[1])
            elif tag == 'END' and current is not None:
                if ('model' in current and len(current['model']) == 11
                        and set(current['momenta']) == set(range(1, 7))
                        and current.get('fks') in range(1, 7)):
                    candidates.append(current)
                current = None
    if not candidates:
        raise ValueError('No complete production momenta/model captures found')
    candidates.sort(key=lambda row: abs(row['net']), reverse=True)
    chosen = candidates[:args.top]
    output = [str(len(chosen))]
    for ordinal, point in enumerate(chosen, 1):
        index = 1 if point['fks'] <= 4 else point['fks']-3
        point.update(probe_id=ordinal, matrix_index=index)
        output += [f'{ordinal} {index}', ' '.join(f'{v:.17e}' for v in point['model'][:5]),
                   ' '.join(f'{v:.17e}' for v in point['model'][5:])]
        output += [' '.join(f'{v:.17e}' for v in point['momenta'][leg]) for leg in range(1, 7)]
    args.output.write_text('\n'.join(output)+'\n')
    metadata = dict(schema='ttw_real_precision_capture_selection_v1', available_points=len(candidates),
        selected=chosen, traces={str(path.resolve()): digest(path) for path in args.traces})
    args.output.with_suffix('.json').write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n')
    print(args.output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    p = sub.add_parser('build')
    p.add_argument('process', type=Path)
    p.add_argument('--subprocess', default=DEFAULT_SUBPROCESS)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--compiler', default='gfortran')
    p.set_defaults(action=build)
    p = sub.add_parser('run')
    p.add_argument('binary', type=Path)
    p.add_argument('input', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.set_defaults(action=run)
    p = sub.add_parser('self-test-input')
    p.add_argument('--output', type=Path, required=True)
    p.set_defaults(action=self_test_input)
    p = sub.add_parser('extract')
    p.add_argument('traces', type=Path, nargs='+')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--top', type=int, default=8)
    p.set_defaults(action=extract)
    args = parser.parse_args()
    args.action(args)


if __name__ == '__main__':
    main()
