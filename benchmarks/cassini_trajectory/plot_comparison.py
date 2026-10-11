"""Plot C++ propagation exports and compute SPK comparison metrics."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LABELS = {"multibody_j2_srp_thermal": "Multibody + complete J2 + rings + 1PN + attitude SRP + thermal recoil"}
COLORS = ["#148c86"]

def make_additional_plots(folder, times, truth, truth_velocity, propagated, velocity, forces, sources, mass):
    days = times/86400
    position_error = (propagated-truth)*1000  # metres
    velocity_error = (velocity-truth_velocity)*1e6  # mm/s
    radial = truth/np.linalg.norm(truth, axis=1)[:,None]
    angular_momentum = np.cross(truth, truth_velocity)
    normal = angular_momentum/np.linalg.norm(angular_momentum, axis=1)[:,None]
    tangent = np.cross(normal, radial)
    basis = np.stack([radial, tangent, normal], axis=1)
    np.testing.assert_allclose(basis @ np.transpose(basis, (0,2,1)),
                               np.broadcast_to(np.eye(3), basis.shape), atol=1e-12)
    rtn = np.einsum("nij,nj->ni", basis, position_error)
    np.testing.assert_allclose(np.linalg.norm(rtn,axis=1), np.linalg.norm(position_error,axis=1),atol=1e-9)
    component_colors = ["#cf6334", "#387bb5", "#8261a8"]
    title = LABELS["multibody_j2_srp_thermal"]

    def save(fig, name):
        fig.savefig(folder / f"{name}.png", dpi=180, bbox_inches="tight")
        fig.savefig(folder / f"{name}.pdf", bbox_inches="tight")
        plt.close(fig)

    def components(ax, values, labels):
        for j, (label, color) in enumerate(zip(labels,component_colors)):
            ax.plot(days, values[:,j], color=color, lw=1.7, label=label)
        ax.axhline(0,color="#52616b",lw=.7,alpha=.5)
        ax.legend(fontsize=9)
        ax.set_xlim(days[0], days[-1])

    fig, axes = plt.subplots(1,3,figsize=(15,5),layout="constrained")
    for ax, (i,j) in zip(axes, [(0,1),(0,2),(1,2)]):
        ax.plot(truth[:,i]/1e6,truth[:,j]/1e6,color="#162f50",lw=3,label="Official SPK")
        ax.plot(propagated[:,i]/1e6,propagated[:,j]/1e6,color=COLORS[0],lw=1.5,ls="--",label="Retained model")
        ax.scatter(truth[0,i]/1e6,truth[0,j]/1e6,color="#162f50",s=25,label="Start")
        ax.scatter(truth[-1,i]/1e6,truth[-1,j]/1e6,color="#162f50",marker="x",s=35,label="End")
        ax.set_xlabel(f"J2000 {'XYZ'[i]} (million km)")
        ax.set_ylabel(f"J2000 {'XYZ'[j]} (million km)")
        ax.set_aspect("equal",adjustable="datalim")
    axes[0].legend(fontsize=8)
    fig.suptitle("Saturn-centered trajectory projections · SPK and retained model overlap",fontsize=14)
    save(fig,"trajectory_projections")

    for values, unit, quantity, filename in [
        (position_error,"m","Position","position_error_components"),
        (velocity_error,"mm/s","Velocity","velocity_error_components")]:
        fig, axes = plt.subplots(2,1,figsize=(11,7),sharex=True,layout="constrained")
        components(axes[0],values,[f"Δ{k}" for k in "XYZ"])
        axes[0].set_ylabel(f"Signed {quantity.lower()} error ({unit})")
        axes[1].plot(days,np.linalg.norm(values,axis=1),color=COLORS[0],lw=2)
        axes[1].set_ylabel(f"Error norm ({unit})");axes[1].set_ylim(bottom=0)
        axes[1].set_xlabel("Elapsed TDB days")
        fig.suptitle(f"{quantity} errors in J2000 · propagated minus SPK\n{title}",fontsize=13)
        save(fig,filename)

    fig, axes = plt.subplots(3,1,figsize=(11,8),sharex=True,layout="constrained")
    for j,(ax,label,color) in enumerate(zip(axes,["Radial","Along-track","Cross-track"],component_colors)):
        ax.plot(days,rtn[:,j],color=color,lw=1.8)
        ax.axhline(0,color="#52616b",lw=.7,alpha=.5)
        ax.set_ylabel(f"{label} error (m)")
    axes[-1].set_xlabel("Elapsed TDB days")
    fig.suptitle("Position divergence in the instantaneous SPK orbital frame",fontsize=14)
    save(fig,"rtn_position_errors")

    fig, axes = plt.subplots(3,1,figsize=(11,8),sharex=True,layout="constrained",
                             gridspec_kw={"height_ratios":[2,2,1]})
    for j,(label,color) in enumerate(zip(["HGA SRP","Thermal recoil"],["#d69822",COLORS[0]])):
        axes[j].plot(days,forces[:,j]*1e9,color=color,lw=1.4)
        axes[j].set_ylabel(f"{label} (nm/s²)")
        axes[j].ticklabel_format(axis="y",style="plain",useOffset=False)
    axes[0].set_ylim(bottom=0)
    codes={"as_flown":0,"reconstructed":1,"short_gap_slerp":2}
    axes[2].step(days,[codes[source] for source in sources],where="mid",color="#546b80",lw=1)
    axes[2].set_yticks([0,1,2],["As-flown","Reconstructed","Short bridge"])
    axes[2].set_ylim(-.3,2.3);axes[2].set_xlabel("Elapsed TDB days")
    fig.suptitle("Individual perturbation magnitudes and attitude sources\nMagnitudes at the exported sampling cadence",fontsize=14)
    save(fig,"perturbation_accelerations")

    fig,ax=plt.subplots(figsize=(8,7),layout="constrained")
    ax.plot(rtn[:,1],rtn[:,0],color="#c5d2db",lw=1,zorder=1)
    dots=ax.scatter(rtn[:,1],rtn[:,0],c=days,cmap="viridis",s=8,zorder=2)
    ax.scatter(rtn[0,1],rtn[0,0],marker="o",color="#162f50",s=45,label="Start",zorder=3)
    ax.scatter(rtn[-1,1],rtn[-1,0],marker="x",color="#cf6334",s=60,label="End",zorder=3)
    ax.set_xlabel("Along-track position error (m)");ax.set_ylabel("Radial position error (m)")
    ax.set_aspect("equal",adjustable="datalim");ax.legend()
    fig.colorbar(dots,ax=ax,label="Elapsed TDB days")
    ax.set_title("Orbital-plane residual path\nComponents in the instantaneous SPK frame")
    save(fig,"orbital_plane_residual")

    fig,axes=plt.subplots(2,2,figsize=(13,8),layout="constrained")
    for ax,values,unit,label in [(axes[0,0],position_error,"m","Position error"),
                                 (axes[0,1],velocity_error,"mm/s","Velocity error")]:
        ax.plot(days,np.linalg.norm(values,axis=1),color=COLORS[0],lw=2)
        ax.set_ylabel(f"{label} norm ({unit})");ax.set_ylim(bottom=0)
    components(axes[1,0],rtn,["Radial","Along-track","Cross-track"])
    axes[1,0].set_ylabel("RTN position error (m)")
    for j,(label,color) in enumerate(zip(["HGA SRP","Thermal recoil"],["#d69822",COLORS[0]])):
        axes[1,1].plot(days,forces[:,j]*1e9,color=color,lw=1.5,label=label)
    axes[1,1].set_ylabel("Acceleration magnitude (nm/s²)")
    axes[1,1].set_ylim(bottom=0);axes[1,1].legend()
    for ax in axes.ravel():ax.set_xlabel("Elapsed TDB days")
    fig.suptitle(f"Cassini · {days[-1]:g}-day SPK comparison · {mass:g} kg\n{title}",fontsize=15)
    save(fig,"diagnostic_overview")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("results/cassini_trajectory"))
    args = parser.parse_args()
    folder = args.input
    metadata = json.loads((folder / "metadata.json").read_text())
    flags = {"indirect_j2":metadata["indirect_j2"], "ring":metadata["rings_enabled"],
             "relativistic":metadata["relativity_enabled"]}
    label = "Multibody + " + ("complete J2" if flags["indirect_j2"] else "direct J2")
    if flags["ring"]: label += " + rings"
    if flags["relativistic"]: label += " + 1PN"
    label += " + attitude SRP + thermal recoil"
    LABELS["multibody_j2_srp_thermal"] = label
    metadata["model_label"] = label
    if not metadata["numerical_validation"]:
        raise ValueError("Study-only output is not a validated retained simulation")
    groups = {}
    with (folder / "trajectory.csv").open() as stream:
        for row in csv.DictReader(stream):
            groups.setdefault(row["model"], []).append(row)
    def array(name, keys):
        return np.array([[float(row[key]) for key in keys] for row in groups[name]])
    times = array("spk", ["elapsed_s"]).ravel()
    assert np.all(np.diff(times) > 0)
    assert times[-1] == metadata["duration_days"] * 86400
    np.testing.assert_allclose(np.diff(times)[:-1], metadata["cadence_s"], rtol=0, atol=1e-9)
    assert 0 < times[-1]-times[-2] <= metadata["cadence_s"]
    for name in groups:
        assert np.array_equal(array(name, ["elapsed_s"]).ravel(), times), name
        assert np.isfinite(array(name, ["x_km", "y_km", "z_km", "position_error_km", "velocity_error_m_s"])).all(), name
    truth = array("spk", ["x_km", "y_km", "z_km"])
    xyz = {name: array(name, ["x_km", "y_km", "z_km"]) for name in LABELS}
    offsets = {name: xyz[name] - truth for name in LABELS}
    for name in LABELS:
        np.testing.assert_allclose(np.linalg.norm(offsets[name], axis=1), array(name, ["position_error_km"]).ravel(), atol=1e-9)
    np.testing.assert_allclose(xyz["multibody_j2_srp_thermal"][0], truth[0], rtol=0, atol=0)
    state_keys = ["x_km", "y_km", "z_km", "vx_km_s", "vy_km_s", "vz_km_s"]
    np.testing.assert_allclose(array("multibody_j2_srp_thermal", state_keys)[0], array("spk", state_keys)[0], rtol=0, atol=0)
    with (folder / "validation.csv").open() as stream:
        validation = list(csv.DictReader(stream))
    np.testing.assert_array_equal([float(row["elapsed_s"]) for row in validation], times)
    refined = np.array([[float(row[key]) for key in ("x_km", "y_km", "z_km")] for row in validation])
    numerical_m = np.linalg.norm(refined-xyz["multibody_j2_srp_thermal"], axis=1)*1000
    if numerical_m.max() > .1:
        raise RuntimeError(f"Integration convergence exceeds 0.1 m: {numerical_m.max()}")
    maneuver_day = metadata["maneuver_elapsed_s"]/86400
    pre = times < metadata["maneuver_elapsed_s"]
    if not pre.any():
        raise ValueError("Arc has no pre-OTM-4 samples; choose an earlier start for this comparison")
    metrics = []
    table_labels = LABELS
    for name, label in table_labels.items():
        error = array(name, ["position_error_km"]).ravel()
        velocity = array(name, ["velocity_error_m_s"]).ravel()
        metrics.append(dict(model=label, initial_km=error[0], pre_burn_rms_km=np.sqrt(np.mean(error[pre]**2)),
                            last_pre_burn_km=error[pre][-1], final_km=error[-1], rms_km=np.sqrt(np.mean(error**2)),
                            max_km=error.max(), final_velocity_m_s=velocity[-1]))
    with (folder / "comparison.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=metrics[0].keys()); writer.writeheader(); writer.writerows(metrics)
    inputs = json.loads(Path("benchmarks/cassini_trajectory/data/gravity_mass_manifest.json").read_text())
    if hashlib.sha256(Path(metadata["gravity_label_path"]).read_bytes()).hexdigest() != inputs["gravity_label_sha256"]:
        raise ValueError("SAT427 label changed; rerun audit_gravity_mass.py")
    for source in inputs["mass_sources"]:
        if hashlib.sha256(Path(source["path"]).read_bytes()).hexdigest() != source["sha256"]:
            raise ValueError("Mass input changed; rerun audit_gravity_mass.py")
    metadata["gravity_mass_inputs"] = inputs
    reference_mass = inputs["mass_kg"]
    metadata["pds_nominal_mass_kg"] = reference_mass
    metadata["mass_offset_from_pds_kg"] = metadata["mass_kg"]-reference_mass
    metadata["inverse_mass_scale_vs_pds"] = reference_mass/metadata["mass_kg"]
    thermal_reference = [metadata[f"thermal_reference_{axis}_km_s2"] for axis in "xyz"]
    nominal_thermal = np.allclose(thermal_reference, [0, 0, -5.45e-12], rtol=0, atol=1e-25)
    recoil_note = ("The reference vector is the paper's later-tour Z mean with X/Y assumed zero, applied as a hypothesis to October 2004."
                  if nominal_thermal else "The reference vector differs from the paper's later-tour Z-only mean; these are diagnostic coefficients, not an October thermal measurement.")
    study_path = Path("benchmarks/cassini_trajectory/data/recoil_study.json")
    if study_path.exists():
        study = json.loads(study_path.read_text())
        selected = study["selected"]
        if (metadata["mass_kg"] == study["mass_kg"] and metadata["gravity_mask"] == study["gravity_mask"]
                and metadata["start_utc"] == study["start_utc"] and metadata["duration_days"] == study["duration_days"]
                and np.allclose(thermal_reference, selected["reference_vector_km_s2"], rtol=0, atol=1e-25)):
            metadata["recoil_study"] = study
            recoil_note += " " + selected["provenance"] + " See [the recoil study](../../benchmarks/cassini_trajectory/RECOIL_UPDATE.md)."
    metadata["max_numerical_difference_m"] = float(numerical_m.max())
    metadata["full_position_rms_km"] = float(metrics[0]["rms_km"])
    with (folder / "force_audit.csv").open() as stream:
        audit = list(csv.DictReader(stream))
    np.testing.assert_array_equal([float(row["elapsed_s"]) for row in audit],times)
    additions = {}
    for key in ["indirect_j2", "ring", "relativistic"]:
        vectors = np.array([[float(row[f"{key}_{axis}_km_s2"]) for axis in "xyz"] for row in audit])
        assert np.isfinite(vectors).all()
        additions[key] = np.linalg.norm(vectors,axis=1)*1e12  # km/s² -> nm/s²
    metadata["added_force_ranges_nm_s2"] = {
        key:dict(minimum=float(value.min()),maximum=float(value.max())) for key,value in additions.items()}
    metadata["attitude_sample_counts"] = {source: sum(row["attitude_source"] == source for row in groups["spk"])
        for source in ("reconstructed", "as_flown", "short_gap_slerp")}
    metadata["attitude_coverage"] = json.loads(Path("benchmarks/cassini_trajectory/attitude_manifest.json").read_text())
    np.testing.assert_allclose(metadata["attitude_reconstructed_seconds"], metadata["attitude_coverage"]["reconstructed_seconds"], rtol=0, atol=1e-6)
    if metadata["attitude_coverage"]["start_utc"] != metadata["start_utc"] or metadata["attitude_coverage"]["end_utc"] != groups["spk"][-1]["utc"]:
        raise ValueError("Run audit_attitude.py for this arc before plotting")
    metadata["last_pre_burn_utc"] = groups["spk"][int(np.where(pre)[0][-1])]["utc"]
    metadata["end_utc"] = groups["spk"][-1]["utc"]
    metadata["kernel_sha256"] = {}
    for kernel in metadata["kernel_paths"]:
        digest = hashlib.sha256()
        with open(kernel, "rb") as stream:
            for chunk in iter(lambda: stream.read(1024*1024), b""): digest.update(chunk)
        metadata["kernel_sha256"][Path(kernel).name] = digest.hexdigest()
    (folder / "summary.json").write_text(json.dumps(metadata, indent=2)+"\n")
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .22, "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(14, 6.5), layout="constrained")
    left = fig.add_subplot(121, projection="3d")
    right = fig.add_subplot(122, projection="3d")
    left.plot(*(truth/1e6).T, color="#162f50", lw=3.5, label="Official reconstructed SPK")
    for (name, label), color in zip(LABELS.items(), COLORS):
        left.plot(*(xyz[name]/1e6).T, color=color, lw=1.4, ls="--", label=label)
        if name != "two_body":
            right.plot(*offsets[name].T, color=color, lw=2, label=label)
    left.scatter(*(truth[0]/1e6), color="#162f50", marker="o", s=40, label="Start")
    left.scatter(*(truth[-1]/1e6), color="#162f50", marker="x", s=45, label="End")
    left.scatter(0, 0, 0, color="#a07b43", s=50, label="Saturn")
    left.set_title("Cassini trajectory relative to Saturn\nSPK and propagated trajectory overlap at this scale")
    right.set_title("Multibody model divergence from SPK\nIntegrated position minus reference position")
    for ax, unit in [(left,"million km"), (right,"km")]:
        ax.set_xlabel(f"J2000 X ({unit})", labelpad=10)
        ax.set_ylabel(f"J2000 Y ({unit})", labelpad=10)
        ax.set_zlabel(f"J2000 Z ({unit})", labelpad=10)
        ax.view_init(elev=25, azim=-55)
        data = np.vstack([truth, *xyz.values(), np.zeros((1,3))])/1e6 if ax is left else np.vstack([v for k,v in offsets.items() if k != "two_body"])
        limits = np.array([data.min(axis=0), data.max(axis=0)])
        center = limits.mean(axis=0); half = max(np.ptp(data, axis=0).max()/2, 1e-6)*1.08
        ax.set_xlim(center[0]-half, center[0]+half); ax.set_ylim(center[1]-half, center[1]+half); ax.set_zlim(center[2]-half, center[2]+half)
        ax.set_box_aspect((1,1,1))
    left.legend(loc="upper left", fontsize=8)
    fig.suptitle(f"Cassini · {metadata['start_utc'][:10]} to {metadata['end_utc'][:10]} · {metadata['duration_days']:g} days · {metadata['mass_kg']:g} kg", fontsize=16)
    fig.savefig(folder / "trajectory_3d.png", dpi=190, bbox_inches="tight")
    fig.savefig(folder / "trajectory_3d.pdf", bbox_inches="tight")
    plt.close(fig)
    # A standalone HTML export lets the user rotate both scientific 3D plots.
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        pass
    else:
        interactive = make_subplots(rows=1, cols=2, specs=[[{"type":"scene"}, {"type":"scene"}]],
                                    subplot_titles=["Saturn-centered trajectory (million km)", "Divergence from SPK (km)"])
        interactive.add_trace(go.Scatter3d(x=truth[:,0]/1e6, y=truth[:,1]/1e6, z=truth[:,2]/1e6,
                              mode="lines", name="Official SPK", line=dict(color="#162f50",width=6)), row=1,col=1)
        for (name,label), color in zip(LABELS.items(),COLORS):
            for col, data in [(1, xyz[name]/1e6),(2,offsets[name])]:
                if col == 2 and name == "two_body":
                    continue
                interactive.add_trace(go.Scatter3d(x=data[:,0], y=data[:,1], z=data[:,2], mode="lines",
                    name=label, legendgroup=name, showlegend=col==1, line=dict(color=color,width=3),
                    customdata=times/86400, hovertemplate="Day %{customdata:.3f}<br>X=%{x:.6f}<br>Y=%{y:.6f}<br>Z=%{z:.6f}<extra>%{fullData.name}</extra>"),row=1,col=col)
        interactive.update_layout(title=f"Cassini · two-week SPK and open-loop dynamics comparison · {metadata['mass_kg']:g} kg",height=700,
            scene=dict(aspectmode="data",xaxis_title="X (million km)",yaxis_title="Y (million km)",zaxis_title="Z (million km)"),
            scene2=dict(aspectmode="data",xaxis_title="ΔX (km)",yaxis_title="ΔY (km)",zaxis_title="ΔZ (km)"),
            legend=dict(orientation="h",y=-.05))
        interactive.write_html(folder / "trajectory_3d.html",include_plotlyjs=True)
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True, layout="constrained")
    for (name,label), color in zip(LABELS.items(),COLORS):
        axes[0].plot(times/86400, array(name,["position_error_km"]).ravel(),color=color,label=label)
        axes[1].plot(times/86400, array(name,["velocity_error_m_s"]).ravel(),color=color,label=label)
    for j, color in enumerate(["#cf6334", "#387bb5", "#8261a8"]):
        axes[2].plot(times/86400,offsets["multibody_j2_srp_thermal"][:,j],color=color,label=f"Δ{'XYZ'[j]}")
    for ax in axes:
        if 0 < maneuver_day < times[-1]/86400:
            ax.axvline(maneuver_day,color="#bd302c",ls=":",label="OTM-4 burn (unmodeled)")
            ax.axvspan(maneuver_day,times[-1]/86400,color="#bd302c",alpha=.06)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("Position error (km)"); axes[1].set_ylabel("Velocity error (m/s)")
    axes[0].set_ylim(bottom=0)
    axes[1].set_ylim(bottom=0)
    axes[2].set_ylabel("Signed position error (km)"); axes[2].set_xlabel("Elapsed TDB days")
    fig.suptitle("Open-loop divergence from Cassini SPK",fontsize=15)
    fig.savefig(folder / "divergence.png",dpi=180); fig.savefig(folder / "divergence.pdf"); plt.close(fig)
    make_additional_plots(folder, times, truth,
        array("spk", ["vx_km_s", "vy_km_s", "vz_km_s"]),
        xyz["multibody_j2_srp_thermal"],
        array("multibody_j2_srp_thermal", ["vx_km_s", "vy_km_s", "vz_km_s"]),
        array("multibody_j2_srp_thermal", ["srp_acceleration_m_s2", "thermal_acceleration_m_s2"]),
        [row["attitude_source"] for row in groups["spk"]],metadata["mass_kg"])
    force_names = {"indirect_j2":"Indirect Saturn J2", "ring":"Ring gravity (direct + origin shape)",
                   "relativistic":"Relative barycentric 1PN"}
    active = [key for key in force_names if flags[key]]
    if active:
        fig, axes = plt.subplots(len(active),1,figsize=(11,3*len(active)),sharex=True,layout="constrained",squeeze=False)
        for key,ax,color in zip(active,axes[:,0],["#8261a8","#cf6334","#387bb5"]):
            ax.plot(times/86400,additions[key],color=color,lw=1.7)
            ax.set_ylabel(f"{force_names[key]}\n(nm/s²)");ax.set_ylim(bottom=0)
        axes[-1,0].set_xlabel("Elapsed TDB days")
        fig.suptitle("Retained gravity additions evaluated along the SPK reference",fontsize=13)
        fig.savefig(folder / "gravity_additions.png",dpi=180)
        fig.savefig(folder / "gravity_additions.pdf");plt.close(fig)
    else:
        for suffix in ["png","pdf"]: (folder / f"gravity_additions.{suffix}").unlink(missing_ok=True)
    burn_inside = 0 <= metadata["maneuver_elapsed_s"] <= times[-1]
    row = metrics[-1]
    lines = ["# Cassini two-week SPK comparison", "",
        f"UTC: {metadata['start_utc']} to {metadata['end_utc']}. Sampling: {metadata['cadence_s']:g} TDB seconds.", "",
        "| Model | Initial position error (km) | Position RMS (km) | Final position error (km) | Maximum position error (km) | Final velocity error (m/s) |",
        "|---|---:|---:|---:|---:|---:|",
        *[f"| {r['model']} | {r['initial_km']:.6f} | {r['rms_km']:.6f} | {r['final_km']:.6f} | {r['max_km']:.6f} | {r['final_velocity_m_s']:.6f} |" for r in metrics], "",
        f"Maximum nominal/refined integration difference: {numerical_m.max():.6f} m.", "",
        "Both position and velocity start from the exact geometric SPK state. No initial offsets, observation updates, or subsequent SPK state resets are applied. Recoil-coefficient calibration, when present, is described below. Only the retained configuration and SPK reference are shown.", "",
        f"SRP uses a one-sided, diffuse HGA flat plate with physical face area A={metadata['area_m2']:g} m², prescribed mass={metadata['mass_kg']:g} kg, nu={metadata['srp_diffuse_nu']:g} and normal-incidence equivalent CR={metadata['cr']:g}. The illuminated projected area is A max(0, n·s), with the HGA normal along spacecraft -Z. Saturn shadowing is retained. Bus, antenna backside and specular reflection are omitted; this is not the complete flight SRP model.", "",
        f"Thermal recoil uses the signed body-frame reference vector [{', '.join(f'{value:.9g}' for value in thermal_reference)}] km/s² at {metadata['thermal_reference_mass_kg']:g} kg, referenced to {metadata['thermal_reference_utc']}. {recoil_note} The reference date within T7 is a convention. The vector scales inversely with prescribed mass and decays with Pu-238 half-life {metadata['thermal_half_life_years']:g} years. It remains active in eclipse.", "",
        f"The NASA PDS nominal October mass is {reference_mass:g} kg from the rounded pre-burn OTM-004 snapshot, supported by surrounding OTM-003/005 entries. This simulation uses {metadata['mass_kg']:g} kg, an offset of {metadata['mass_offset_from_pds_kg']:+g} kg from that source. A nonzero offset is a requested sensitivity scenario, not a revised mass measurement. Both non-gravitational terms scale by {metadata['inverse_mass_scale_vs_pds']:.8f} relative to the nominal mass at fixed geometry. The table has roughly 10 kg reporting resolution and does not resolve continuous consumption. Huygens was still attached. [NASA PDS mass table](https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/Cassini/files/Cassini%20Spacecraft%20Mass%20Table.csv). Optical coefficients remain fixed. [Di Benedetto, Iess & Roth (2009), §2 and Appendix](https://issfd.org/ISSFD_2009/OrbitDeterminationI/DiBENEDETTO.pdf).", "",
        f"Attitude coverage: {metadata['attitude_coverage']['reconstructed_seconds']/86400:.5f} days reconstructed, {metadata['attitude_coverage']['as_flown_seconds']/86400:.5f} days from predicted CKs generated from the as-flown sequence. A {metadata['attitude_coverage']['gap_seconds']:.6f} s gap is bridged by SLERP, with an enforced 10 s limit and no extrapolation. Reconstructed kernels have load priority. See the tracked attitude_manifest.json for URLs, hashes and gap times.", "",
        ("OTM-4 lies inside this arc and is unmodeled." if burn_inside else "This arc ends before OTM-4: the mass table lists 23 October 06:16 UTC, and the weekly event report gives 07:29:28 UTC. Both are later than the endpoint; the earlier time is used for pre-burn checks. Absence of all small attitude-control firings has not been established."),
        "[NASA OTM-4 event report](https://science.nasa.gov/missions/cassini/significant-event-report-for-week-ending-10292004/).", "",
        f"Multibody gravity includes Sun, eight other planetary-system barycenters and ten Saturn moons, with indirect origin acceleration. SAT427.23’s paired SPK label supplies Saturn/moon GMs, Jupiter/Uranus/Neptune system GMs, J2={metadata['j2']:.16g}, gravity reference radius {metadata['saturn_radius_km']:g} km and the three-term analytic pole series (TDB Julian centuries; RA sine, declination cosine). The physical shadow radius is separately {metadata['saturn_shadow_radius_km']:g} km. SAT427’s solar GM includes inner planets; explicit barycenter GMs 1–4 are subtracted to avoid double counting. Remaining planetary GMs use the generic DE440 kernel. [SAT427 SPK release label](https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/spk/200128RU_SCPSE_04251_04303.bsp.lbl), [pole-series convention](https://doi.org/10.3847/1538-3881/ac90c9). Higher Saturn harmonics, moon shadows, bus/backside SRP and maneuvers are omitted. Numerical convergence does not validate these physical assumptions.", "",
        f"Gravity mask: {metadata['gravity_mask']}; selected default: {metadata['retained_default_gravity_mask']}. Indirect J2 enabled: {flags['indirect_j2']}; rings enabled: {flags['ring']}; relativity enabled: {flags['relativistic']}. The default minimized sampled whole-arc position RMS among eight combinations tested at 3220 kg with the original nominal Z-only recoil and all other settings fixed. Changes to mass or recoil keep this gravity selection without repeating its ranking. Other masks are diagnostic overrides. This SPK-arc selection is not independent physical validation; disabling a required origin correction does not make it physically unnecessary. See [the ablation audit](../../benchmarks/cassini_trajectory/GRAVITY_UPDATE.md).", "",
        "The available J2 origin correction sums the moons' prescribed J2 accelerations with positive GM weighting. The available relativity model uses Moyer Eq. 4-26 with beta=gamma=1, 20 point-mass sources and barycentric velocities, subtracting Saturn's 1PN acceleration from Cassini's. [JPL formulation, §§4.4.1 and 4.4.6](https://descanso.jpl.nasa.gov/monograph/series2/Descanso2_all.pdf).", "",
        "When enabled, the separate C/B/A rings use SAT427 GM values 0.039/0.784/0.300 km³/s² and its listed inner/outer radii, uniformly distributed over thin annuli with the same pole. Their exterior Newton fields are evaluated through degree 80. The direct monopole and shape-driven origin reaction are included with a composite central mass. Uniform density and concentric co-acceleration are approximations, not a verified reconstruction of SAT427's ring implementation. Relativistic ring/multipole and spin terms are omitted. [Disk/annulus Newtonian field reference](https://ntrs.nasa.gov/citations/19820045807).", "",
        "| Assessed gravity term on SPK path | Enabled in retained model | Minimum (nm/s²) | Maximum (nm/s²) |",
        "|---|---|---:|---:|",
        *[f"| {force_names[key]} | {'Yes' if flags[key] else 'No'} | {values.min():.6f} | {values.max():.6f} |" for key,values in additions.items()], "",
        "## Plot gallery", "",
        "[Rotatable 3D trajectory and divergence](trajectory_3d.html)", "",
        "![Diagnostic overview](diagnostic_overview.png)", "",
        "![3D trajectory](trajectory_3d.png)", "", "![Divergence](divergence.png)", "",
        "![Trajectory projections](trajectory_projections.png)", "",
        "![Position errors](position_error_components.png)", "",
        "![Velocity errors](velocity_error_components.png)", "",
        "![Radial, along-track and cross-track position errors](rtn_position_errors.png)", "",
        "RTN axes are defined from the SPK position and angular momentum at each epoch. Components are propagated position minus SPK position projected into that instantaneous frame.", "",
        "![SRP, thermal recoil and attitude sources](perturbation_accelerations.png)", "",
        *( ["![Retained gravity additions](gravity_additions.png)", ""] if active else [] ),
        "Acceleration curves show individual force magnitudes at the output cadence. The attitude strip labels output samples; exact coverage durations and the short bridge are recorded above.", "",
        "![Residual path in the orbital plane](orbital_plane_residual.png)", "",
        "Each figure is also available as a PDF next to its PNG.", ""]
    report = "\n".join(lines)
    (folder / "COMPARISON.md").write_text(report)
    print(report.split("![3D")[0])

if __name__ == "__main__":
    main()
