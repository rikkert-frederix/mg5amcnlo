"""Analytic figure-value checks and an explicitly synthetic layout fixture."""
import copy
from pathlib import Path
import tempfile
import unittest

import numpy as np

from campaign import digest
from main_figures import BANDS, band_curve, row_values, run, series
from pilot_report import scales


def fixture_row(value,error=None):
    points={'central':value}
    points.update({point:value*(1+.08*np.log(point[0])+.03*np.log(point[1])+
                              .02*np.log(point[2])+.01*np.log(point[3])) for point in scales.POINTS})
    return dict(scale=scales.describe(points),nominal_retraining_mc_error=abs(value)*.015 if error is None else error)


class SyntheticFigures:
    banner='SYNTHETIC LAYOUT TEST DATA — no physics predictions'
    scope='synthetic_plot_layout_only'
    inputs={str(Path(__file__).resolve()):digest(Path(__file__))}

    def rates(self,label,charge='plus',config='R04_b25',derived=False):
        base={'fiducial_2b':.1,'2b_0extra':.07,'2b_1extra':.025,'2b_2extra':.004,
              '2b_3plus_extra':.00005,'two_b_given_one_b':.72,'veto_0extra':.7}[label]
        if charge=='minus':base*=.6
        base*=dict(R04_b25=1.,R04_b30=.97,R04_b40=.9,R03_b25=.95,R05_b25=1.04)[config]
        rows={variant:fixture_row(base*factor) for variant,factor in (('S',1.),('Pi',1.025),('P',1.15),('PiD',.96))}
        if label in ('2b_2extra','2b_3plus_extra'):rows['S']=fixture_row(0.)
        if label=='2b_3plus_extra':rows['PiD']=fixture_row(0.)
        return rows

    def histogram(self,local,charge='plus',charge_comparison=False):
        high={3:800,4:600,5:5,6:np.pi,7:6,8:1600,11:400,12:400,16:400}[local]
        boundaries=np.linspace(0,high,9)
        edges=np.column_stack([boundaries[:-1],boundaries[1:]])
        fractions=np.exp(-np.arange(8)/2.)
        fractions=.85*fractions/fractions.sum()
        rows=[]
        for i,fraction in enumerate(fractions):
            record={variant:fixture_row(fraction*(1+tilt*(i-2))) for variant,tilt in
                    (('S',0.),('Pi',.018),('P',-.025),('PiD',.009))}
            record['Pi_over_S']=fixture_row(1+.018*(i-2))
            rows.append(record)
        detail=dict(title='synthetic observable '+str(local),edges=edges.tolist(),
                    normalization_rate_bins=['fiducial_2b' if local in (11,12,16) else 'fiducial_1b'])
        if charge_comparison:
            detail['normalized']=[{variant:{'plus_shape_over_minus_shape':fixture_row(1+factor*(i-2))}
                for variant,factor in (('S',.01),('Pi',.015))} for i in range(8)]
        else:
            detail['normalized']=rows
            detail['absolute']=copy.deepcopy(rows)
            for row in detail['absolute']:
                for variant in ('S','Pi','P','PiD'):
                    row[variant]=fixture_row(row[variant]['scale']['central']*.1)
                if local==16:row['S']=fixture_row(0.)
        return detail


class TestMainFigures(unittest.TestCase):
    def test_band_points_are_verified_before_plotting(self):
        row=fixture_row(10.,.2)
        for band in BANDS:
            values=row_values(row,band)
            self.assertEqual(values[0],10.)
            self.assertEqual(values[1],.2)
            self.assertLessEqual(values[2],10.)
            self.assertGreaterEqual(values[3],10.)
        bad=copy.deepcopy(row)
        bad['scale']['combined63']['envelope'][0]=0.
        with self.assertRaisesRegex(ValueError,'declared matching scale points'):
            row_values(bad)
        bad=copy.deepcopy(row)
        del bad['scale']['points']['2,0.5,1,1']
        with self.assertRaisesRegex(ValueError,'complete signed scale grid'):
            row_values(bad)

    def test_density_conversion_preserves_bin_fractions_and_mc_errors(self):
        detail=dict(edges=[[0.,2.],[2.,5.]],normalized=[{'S':fixture_row(.2,.02)},{'S':fixture_row(.3,.03)}])
        edges,value=series(detail,'S')
        np.testing.assert_allclose(value[:,:2],[[.1,.01],[.1,.01]])
        self.assertAlmostEqual(np.sum(value[:,0]*np.diff(edges,axis=1)[:,0]),.5)
        unchanged=series(detail,'S',density=False)[1]
        np.testing.assert_allclose(unchanged[:,0],[.2,.3])

    def test_unfinished_data_produces_no_figure_directory(self):
        import json
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'pending.json'
            path.write_text(json.dumps(dict(status='prepared; technical/reference release required')))
            target=Path(folder)/'figures'
            with self.assertRaisesRegex(ValueError,'complete full-flavour main reduction'):
                run(path,target)
            self.assertFalse(target.exists())

    def test_ratio_step_does_not_add_zero_at_range_boundaries(self):
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots()
        band_curve(ax,np.array([[0.,1.],[1.,2.]]),np.array([[1.,.01,.95,1.05],[1.1,.01,1.05,1.15]]),
                   'blue','ratio')
        self.assertGreaterEqual(np.min(ax.lines[0].get_ydata()),1.)
        self.assertGreater(ax.get_ylim()[0],.8)
        plt.close(fig)


if __name__=='__main__':
    unittest.main()
