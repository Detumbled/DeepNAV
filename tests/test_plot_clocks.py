"""Clock CSV/plotter contract checks; run: python3 -m unittest discover -s tests -p test_plot_clocks.py."""

import importlib.util
import contextlib
import io
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import matplotlib

# Keep all plotting tests headless, including tests that mock --show.
matplotlib.use("Agg")


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("plot_clocks", ROOT / "plot_clocks.py")
plotter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(plotter)


class ClockPlots(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def csv(self, content, name="history.csv"):
        path = self.directory / name
        path.write_text(content)
        return path

    def test_invalid_csv(self):
        header = "time_s,bias_s,fractional_frequency\n"
        invalid = [
            "time_s,bias_s\n0,0\n1,0\n",
            "time_s,bias_s,bias_s,fractional_frequency\n0,0,0,0\n1,0,0,0\n",
            header,
            header + "0,0,0\n",
            header + "0,0,0\n0,1,1\n",
            header + "1,0,0\n0,1,1\n",
            header + "-1,0,0\n1,1,1\n",
            header + "0,0,0\n1,nan,1\n",
            header + "0,0,0\n1,inf,1\n",
            header + "0,0,0\n1,,1\n",
            header + "0,0,0\n1,word,1\n",
            header + "0,0,0\n\n1,0,0\n",
            header + "0,0,0,7\n1,0,0,7\n",
            "time_s,bias_s,fractional_frequency,sigma_bias_s\n0,0,0,0\n1,0,0,0\n",
            "time_s,bias_s,fractional_frequency,sigma_bias_s,sigma_fractional_frequency\n"
            "0,0,0,0,0\n1,0,0,-1,0\n",
            "time_s,bias_s,fractional_frequency,sigma_bias_s,sigma_fractional_frequency\n"
            "0,0,0,0,0\n1,0,0,0,-1\n",
        ]
        for content in invalid:
            with self.subTest(content=content), self.assertRaises(ValueError):
                plotter.load_clock_csv(self.csv(content))
        with self.assertRaises(ValueError):
            plotter.load_clock_csv(self.directory / "missing.csv")

    def test_precision_and_selection(self):
        bias = np.nextafter(1e-9, 2e-9)
        path = self.csv("time_s,bias_s,fractional_frequency\n"
                        f"0,{bias:.17g},-1e-12\n2,2e-9,3e-12\n")
        frame = plotter.load_clock_csv(path)
        self.assertEqual(frame["bias_s"].iloc[0], bias)
        self.assertTrue(all(str(dtype) == "float64" for dtype in frame.dtypes))
        full = plotter.load_clock_csv(self.csv("time_s,bias_s,fractional_frequency\n" +
            "".join(f"{i},0,0\n" for i in range(100))))
        shown = plotter.display_frame(full, 7)
        self.assertEqual(len(shown), 7)
        self.assertEqual(shown.iloc[0]["time_s"], 0)
        self.assertEqual(shown.iloc[-1]["time_s"], 99)
        self.assertEqual(len(full), 100)
        self.assertIs(plotter.display_frame(full, 0), full)

    def test_plot_units_thresholds_and_formats(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        frame = plotter.load_clock_csv(self.csv(
            "time_s,bias_s,fractional_frequency,sigma_bias_s,sigma_fractional_frequency\n"
            "0,1e-9,1e-12,2e-9,1e-15\n86400,-2e-9,2e-12,3e-9,2e-15\n"))
        second = plotter.load_clock_csv(self.csv(
            "time_s,bias_s,fractional_frequency,sigma_bias_s,sigma_fractional_frequency\n"
            "0,0,0,1e-11,1e-16\n43200,1e-9,1e-15,2e-11,2e-16\n", "second.csv"))
        mean = frame[list(plotter.REQUIRED)].copy()
        mean["bias_s"] = [5e-10, -8e-10] # Explicit expected bias, different from truth.
        second_mean = second[list(plotter.REQUIRED)].copy()
        original_subplots = plt.subplots
        captured = []

        def capture(*args, **kwargs):
            result = original_subplots(*args, **kwargs)
            captured.append(result)
            return result

        with patch.object(plt, "subplots", side_effect=capture):
            for extension in ("png", "pdf", "svg"):
                path = self.directory / f"plot.{extension}"
                outputs = plotter.plot_histories([frame, second], ["A", "B"], path,
                    time_unit="days", k=4, threshold_m=1, mean_frames=[mean, second_mean])
                self.assertEqual(len(outputs), 3)
                for output in outputs:
                    self.assertGreater(output.stat().st_size, 100)
        for offset in range(0, len(captured), 3):
            fig, axes = captured[offset]
            sigma_fig, sigma_axes = captured[offset + 1]
            diagnostic_fig, diagnostic_axes = captured[offset + 2]
            for figure in (fig, sigma_fig, diagnostic_fig):
                self.assertFalse(plt.fignum_exists(figure.number))
            np.testing.assert_array_equal(axes[0, 0].lines[0].get_xdata(), [0, 1])
            np.testing.assert_array_equal(axes[0, 1].lines[0].get_xdata(), [0, 0.5])
            np.testing.assert_allclose(diagnostic_axes[0, 0].lines[0].get_ydata(), [1, -2])
            np.testing.assert_allclose(axes[0, 0].lines[0].get_ydata(), plotter.C * np.array([1e-9, -2e-9]))
            np.testing.assert_allclose(sigma_axes[0, 0].lines[0].get_ydata(), 4 * plotter.C * np.array([2e-9, 3e-9]))
            np.testing.assert_allclose(axes[1, 0].lines[0].get_ydata(),
                plotter.C * (np.array([5e-10, 8e-10]) + 4 * np.array([2e-9, 3e-9])))
            np.testing.assert_array_equal(axes[0, 0].lines[-2].get_ydata(), [1, 1])
            np.testing.assert_array_equal(axes[0, 0].lines[-1].get_ydata(), [-1, -1])
            self.assertNotEqual(axes[0, 0].get_ylim(), axes[0, 1].get_ylim())
            self.assertEqual(sigma_axes[0, 0].get_ylim(), sigma_axes[0, 1].get_ylim())
            self.assertEqual(axes[0, 1].get_xlim(), (0, 1))

    def test_contact_before_display_selection(self):
        frame = plotter.load_clock_csv(self.csv("time_s,bias_s,fractional_frequency\n" +
            "".join(f"{i},{2 if i == 1 else 0},0\n" for i in range(100))))
        shown = plotter.display_frame(frame, 2)
        self.assertTrue((shown["bias_s"] == 0).all())
        self.assertEqual(plotter.first_recorded_contact(frame, frame["bias_s"], 1), (1, 1, 0))
        self.assertIsNone(plotter.first_recorded_contact(frame, frame["bias_s"], 3))
        self.assertEqual(plotter.first_recorded_contact(frame, -frame["bias_s"], 2), (1, 1, 0))

    def test_missing_or_mismatched_means(self):
        frame = plotter.load_clock_csv(self.csv("time_s,bias_s,fractional_frequency\n0,0,0\n1,0,0\n"))
        with self.assertRaisesRegex(ValueError, "requires covariance"):
            plotter.plot_histories([frame], ["A"], self.directory / "bad.png", mean_frames=[frame])
        covariance = frame.copy()
        covariance["sigma_bias_s"] = 0
        covariance["sigma_fractional_frequency"] = 0
        mismatch = frame.copy()
        mismatch["time_s"] = [0, 2]
        with self.assertRaisesRegex(ValueError, "timestamps must match"):
            plotter.plot_histories([covariance], ["A"], self.directory / "bad.png", mean_frames=[mismatch])
        outputs = plotter.plot_histories([frame], ["A"], self.directory / "no_mean.png")
        self.assertEqual(len(outputs), 3) # No mean or covariance is invented.

    def test_cli_validation(self):
        for options in (["--max-points", "1"], ["--max-points", "-1"], ["--k", "nan"],
                        ["--k", "-1"], ["--threshold-m", "0"], ["--threshold-m", "inf"],
                        ["--labels", "one", "two"], ["--output", "plot.jpg"]):
            with self.subTest(options=options), contextlib.redirect_stderr(io.StringIO()), \
                    self.assertRaises(SystemExit) as raised:
                plotter.main(["missing.csv"] + options)
            self.assertEqual(raised.exception.code, 2)

    def test_no_arguments_uses_setup(self):
        sentinel = object()
        with patch.object(plotter, "load_clock_csv", return_value=sentinel) as load, \
                patch.object(plotter, "plot_histories", return_value=[self.directory / "plot.png"]) as draw, \
                contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(plotter.main([]), 0)
        directory = plotter.SETUP["output_directory"]
        self.assertEqual([call.args[0] for call in load.call_args_list],
                         [directory / name for key in ("csv_files", "mean_csv", "allan_csv")
                          for name in plotter.SETUP[key]])
        args = draw.call_args.args
        self.assertEqual(args[1], plotter.SETUP["labels"])
        self.assertEqual(args[2], directory / plotter.SETUP["output"])
        self.assertEqual(args[3], "days")
        self.assertEqual(args[4:8], (3.0, 1.0, 10000, True))
        with patch.dict(plotter.SETUP, {"k": float("nan")}), \
                contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            plotter.main([])

    def test_setup_optional_settings(self):
        with patch.dict(plotter.SETUP, {"mean_csv": None, "allan_csv": None,
                                       "threshold_m": None, "show": False}):
            args = plotter.setup_arguments()
        for option in ("--mean-csv", "--allan-csv", "--threshold-m", "--show"):
            self.assertNotIn(option, args)

    def test_show_keeps_all_windows_open_until_show_returns(self):
        import matplotlib.pyplot as plt
        frame = plotter.load_clock_csv(self.csv("time_s,bias_s,fractional_frequency\n0,0,0\n1,0,0\n"))
        with patch.object(plt, "show", side_effect=lambda: self.assertEqual(len(plt.get_fignums()), 3)) as show:
            plotter.plot_histories([frame], ["A"], self.directory / "windows.png", show=True)
            show.assert_called_once()
        self.assertFalse(plt.get_fignums())

    def test_allan_curve_noise_drift_and_sampling(self):
        import pandas as pd
        for qb, qy in ((1e-24, 0), (0, 3e-26), (1e-24, 3e-26)):
            rng = np.random.default_rng(123)
            z0, z1, z2 = rng.standard_normal((3, 131072))
            dy = np.sqrt(qy) * z1
            db = np.sqrt(qb) * z0 + np.sqrt(qy) * (z1 / 2 + z2 / np.sqrt(12))
            y = np.r_[0, np.cumsum(dy)]
            bias = np.r_[0, np.cumsum(y[:-1] + db)]
            frame = pd.DataFrame({"time_s": np.arange(len(bias), dtype=float), "bias_s": bias})
            tau, adev = plotter.overlapping_allan_curve(frame)
            short = tau <= 256
            theoretical = np.sqrt(qb / tau + qy * tau / 3)
            np.testing.assert_allclose(adev[short] / theoretical[short], 1, rtol=0.20)
            self.assertLessEqual(tau[-1], (len(bias) - 1) / 10)
        t = np.arange(1001, dtype=float) * 10
        drift = 2e-15
        frame = pd.DataFrame({"time_s": t, "bias_s": 0.5 * drift * t*t})
        tau, adev = plotter.overlapping_allan_curve(frame)
        np.testing.assert_allclose(adev, drift * tau / np.sqrt(2), rtol=1e-9)
        frame["bias_s"] = 1e-8 + 2e-11*t
        self.assertLess(plotter.overlapping_allan_curve(frame)[1].max(), 1e-22)
        frame.loc[1000, "time_s"] += 1
        with self.assertRaisesRegex(ValueError, "uniformly sampled"):
            plotter.overlapping_allan_curve(frame)
        with self.assertRaisesRegex(ValueError, "three samples"):
            plotter.overlapping_allan_curve(frame.iloc[:2])

    def test_allan_log_axes_and_four_windows(self):
        import pandas as pd
        import matplotlib.pyplot as plt
        t = np.arange(101, dtype=float)
        frames = [pd.DataFrame({"time_s": t, "bias_s": d*t*t,
                                "fractional_frequency": 2*d*t}) for d in (1e-15, 1e-17)]
        original = plt.subplots
        captured = []
        def capture(*args, **kwargs):
            result = original(*args, **kwargs)
            captured.append(result)
            return result
        with patch.object(plt, "subplots", side_effect=capture), \
                patch.object(plt, "show", side_effect=lambda: self.assertEqual(len(plt.get_fignums()), 4)):
            paths = plotter.plot_histories(frames, ["Local", "DSAC"], self.directory / "four.png",
                                           max_points=2, show=True, allan_frames=frames, q_bias_s=[0, 0],
                                           drift_per_s=[2e-15, 2e-17], uso_inset_index=0)
        self.assertEqual(len(paths), 4)
        self.assertEqual(paths[-1].name, "four_allan.png")
        _, axes = captured[-1]
        self.assertEqual(axes.shape, (1, 2))
        for column in range(2):
            self.assertEqual(axes[0, column].get_xscale(), "log")
            self.assertEqual(axes[0, column].get_yscale(), "log")
            tau, adev = plotter.overlapping_allan_curve(frames[column])
            np.testing.assert_array_equal(axes[0, column].lines[0].get_xdata(), tau)
            np.testing.assert_allclose(axes[0, column].lines[0].get_ydata(), adev)
            np.testing.assert_allclose(axes[0, column].lines[1].get_ydata(), adev)
            self.assertEqual(axes[0, column].lines[1].get_linestyle(), "--")
        self.assertEqual(len(captured[0][1][0, 0].child_axes), 1)
        self.assertFalse(plt.get_fignums())

    def test_white_fm_aging_theory_and_budget(self):
        import pandas as pd
        t = np.arange(131073, dtype=float)
        for qb, drift in ((2.5e-25, 1e-10/86400), (7.776e-25, 3e-16/86400)):
            rng = np.random.default_rng(123)
            bias = np.r_[0, np.cumsum(np.sqrt(qb)*rng.standard_normal(len(t)-1))] + 0.5*drift*t*t
            tau, adev = plotter.overlapping_allan_curve(pd.DataFrame({"time_s": t, "bias_s": bias}))
            short = tau <= 256
            np.testing.assert_allclose(adev[short] / plotter.white_fm_adev(tau[short], qb, drift), 1, rtol=0.2)
        crossing = plotter.analytical_budget_crossing(1, 7.776e-25, 3e-16/86400) / 86400
        self.assertAlmostEqual(crossing, 8.87027, places=5)
        stochastic = plotter.analytical_budget_crossing(1, 7.776e-25, 0) / 86400
        self.assertGreater(stochastic, crossing)
        self.assertIsNone(plotter.analytical_budget_crossing(1, 0, 0))
        self.assertAlmostEqual(plotter.analytical_budget_crossing(1, 0, 1e-10/86400), 2400.83, places=2)

    def test_cpp_csv_to_python(self):
        binary = ROOT / "build-clang/test_clock_history"
        if not binary.exists():
            self.skipTest("Build test_clock_history to run the C++/Python integration check")
        subprocess.run([str(binary), str(self.directory)], check=True, capture_output=True, text=True)
        local = plotter.load_clock_csv(self.directory / "local.csv")
        dsac = plotter.load_clock_csv(self.directory / "dsac.csv")
        self.assertEqual(len(local), len(dsac))
        full_local = plotter.load_clock_csv(self.directory / "local_allan.csv")
        full_dsac = plotter.load_clock_csv(self.directory / "dsac_allan.csv")
        np.testing.assert_array_equal(full_local["time_s"], full_dsac["time_s"])
        self.assertEqual(local["time_s"].iloc[-1], full_local["time_s"].iloc[-1])
        self.assertEqual(dsac["time_s"].iloc[-1], full_dsac["time_s"].iloc[-1])
        indices = np.searchsorted(full_local["time_s"], local["time_s"])
        np.testing.assert_array_equal(full_local.iloc[indices].to_numpy(), local.to_numpy())
        plotter.overlapping_allan_curve(full_local)
        self.assertIn("sigma_bias_s", local)
        self.assertIn("sigma_bias_s", dsac)
        np.testing.assert_array_equal(local["time_s"], dsac["time_s"])
        local_mean = plotter.load_clock_csv(self.directory / "local_mean.csv")
        t = local_mean["time_s"]
        np.testing.assert_allclose(local_mean["bias_s"], 0.5*(1e-10/86400)*t*t)
        np.testing.assert_allclose(local["sigma_bias_s"], np.sqrt(2.5e-25 * local["time_s"]), rtol=1e-10)
        np.testing.assert_allclose(dsac["sigma_bias_s"], np.sqrt(7.776e-25 * dsac["time_s"]), rtol=1e-10)
        np.testing.assert_allclose(dsac["fractional_frequency"], (3e-16 / 86400) * dsac["time_s"], rtol=1e-10)
        dsac_mean = plotter.load_clock_csv(self.directory / "dsac_mean.csv")
        np.testing.assert_allclose(dsac_mean["bias_s"], 0.5 * (3e-16 / 86400) * t*t, rtol=1e-10)
        result = plotter.main([str(self.directory / "local.csv"), str(self.directory / "dsac.csv"),
                               "--labels", "USO (aging + white FM)", "DSAC (day-matched white FM)",
                               "--mean-csv", str(self.directory / "local_mean.csv"), str(self.directory / "dsac_mean.csv"),
                               "--max-points", "12", "--output", str(self.directory / "comparison.png")])
        self.assertEqual(result, 0)
        self.assertTrue((self.directory / "comparison_allan.png").is_file())


if __name__ == "__main__":
    unittest.main()
