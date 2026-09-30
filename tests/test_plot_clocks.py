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

    def test_show_keeps_all_windows_open_until_show_returns(self):
        import matplotlib.pyplot as plt
        frame = plotter.load_clock_csv(self.csv("time_s,bias_s,fractional_frequency\n0,0,0\n1,0,0\n"))
        with patch.object(plt, "show", side_effect=lambda: self.assertEqual(len(plt.get_fignums()), 3)) as show:
            plotter.plot_histories([frame], ["A"], self.directory / "windows.png", show=True)
            show.assert_called_once()
        self.assertFalse(plt.get_fignums())

    def test_cpp_csv_to_python(self):
        binary = ROOT / "build-clang/test_clock_history"
        if not binary.exists():
            self.skipTest("Build test_clock_history to run the C++/Python integration check")
        subprocess.run([str(binary), str(self.directory)], check=True, capture_output=True, text=True)
        local = plotter.load_clock_csv(self.directory / "local.csv")
        dsac = plotter.load_clock_csv(self.directory / "dsac.csv")
        self.assertEqual(len(local), 61)
        self.assertEqual(len(dsac), 61)
        self.assertEqual(local["time_s"].iloc[-1], 86400)
        self.assertEqual(dsac["time_s"].iloc[-1], 86400)
        self.assertIn("sigma_bias_s", local)
        self.assertIn("sigma_bias_s", dsac)
        np.testing.assert_array_equal(local["time_s"], dsac["time_s"])
        local_mean = plotter.load_clock_csv(self.directory / "local_mean.csv")
        t = local_mean["time_s"]
        np.testing.assert_allclose(local_mean["bias_s"], 2e-9 + 1e-11*t + 0.5*2e-15*t*t)
        np.testing.assert_allclose(local["sigma_bias_s"].iloc[-1],
                                   np.sqrt(1e-24 * 86400 + 3e-26 * 86400**3 / 3))
        result = plotter.main([str(self.directory / "local.csv"), str(self.directory / "dsac.csv"),
                               "--labels", "Local (synthetic)", "DSAC (custom)",
                               "--mean-csv", str(self.directory / "local_mean.csv"), str(self.directory / "dsac_mean.csv"),
                               "--max-points", "12", "--output", str(self.directory / "comparison.png")])
        self.assertEqual(result, 0)


if __name__ == "__main__":
    unittest.main()
