import json
from pathlib import Path
import tempfile
import unittest

from campaign import digest
from main_tables import VARIANTS, band_tex, estimate_tex, render, run, texts
from test_main_figures import fixture_row


class SyntheticTables:
    scope='synthetic_plot_layout_only'
    inputs={str(Path(__file__).resolve()):digest(Path(__file__))}

    def __init__(self):
        self.report=dict(charge_rates={self.rate_title().replace('W+','W+/W-'):{'fiducial_2b':{
            variant:dict(plus_over_minus=fixture_row(1/.6,.03),asymmetry=fixture_row(.25,.009))
            for variant in VARIANTS}}})

    def rate_title(self):return 'R04_b25 W+ rates: synthetic'

    def rates(self,selection,charge='plus',config='R04_b25',derived=False):
        base=dict(fiducial_1b=.0014,fiducial_2b=.001,acceptance_1b=.5,acceptance_2b=.36)[selection]
        if charge=='minus':base*=.95 if derived else .6
        factors=dict(LO=.9,P=1.15,D=.95,S=1.,PiD=.975,Pi=1.025)
        rows={variant:fixture_row(base*factors[variant]) for variant in VARIANTS}
        rows['Pi_minus_S']=fixture_row(base*.025,base*.025)
        rows['Pi_over_S']=fixture_row(1.025,.025)
        rows['Pi_over_S']['pdf']=dict(error_symmetric=.003)
        return rows


class TestMainTables(unittest.TestCase):
    def test_rate_units_and_matched_error_precision(self):
        self.assertEqual(estimate_tex(fixture_row(.001,.000015),1000.),r'$1.000\pm0.015$')
        self.assertEqual(estimate_tex(fixture_row(0.,0.)),r'$0\pm0$')
        row=fixture_row(1.)
        row['nominal_retraining_mc_error']=None
        self.assertEqual(estimate_tex(row),r'$\mathrm{undefined}$')
        row['nominal_retraining_mc_error']=-.1
        with self.assertRaisesRegex(ValueError,'Negative MC error'):
            estimate_tex(row)
        self.assertEqual(band_tex(fixture_row(0.),'combined63',relative=True),r'$\mathrm{undefined}$')

    def test_complete_table_content_and_no_overwrite(self):
        data=SyntheticTables()
        content=texts(data)
        self.assertEqual(set(content),{'main_rates.tex','main_responses.tex'})
        self.assertIn('SYNTHETIC LAYOUT TEST',content['main_rates.tex'])
        self.assertIn(r'$1.000\pm0.015$',content['main_rates.tex'])
        self.assertIn(r'$\delta_{\rm PDF}R_\Pi$',content['main_responses.tex'])
        self.assertIn(r'\rho_X',content['main_responses.tex'])
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'tables'
            render(data,output)
            record=json.loads((output/'manifest.json').read_text())
            self.assertIn('synthetic',record['status'])
            for path,sha in record['outputs'].items():self.assertEqual(digest(path),sha)
            with self.assertRaisesRegex(ValueError,'overwrite'):
                render(data,output)

    def test_pending_main_is_rejected_before_writing(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'pending.json'
            path.write_text(json.dumps(dict(status='prepared')))
            output=Path(directory)/'tables'
            with self.assertRaisesRegex(ValueError,'complete full-flavour main reduction'):
                run(path,output)
            self.assertFalse(output.exists())


if __name__=='__main__':
    unittest.main()
