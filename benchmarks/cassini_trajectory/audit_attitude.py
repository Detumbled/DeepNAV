"""Audit Cassini CK provenance and exact interval coverage for the benchmark arc."""
import argparse
import hashlib
import json
from pathlib import Path
import spiceypy as s

FILES = ["04256_04292py_as_flown.bc", "04292_04320py_as_flown.bc",
         "04281_04286ra.bc", "04286_04291ra.bc", "04291_04296ra.bc",
         "04296_04301ra.bc", "04293_04298ra.bc"]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-utc", default="2004-10-09T00:00:00")
    parser.add_argument("--days", type=float, default=14)
    parser.add_argument("--output", type=Path, default=Path("benchmarks/cassini_trajectory/attitude_manifest.json"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    for name in ["naif0012.tls", "Cassini/sclk/cas00172.tsc", "Cassini/fk/cas_v40.tf"]:
        s.furnsh(str(root / "Kernels" / name))
    start=s.str2et(args.start_utc); end=start+args.days*86400
    all_windows=s.utils.support_types.SPICEDOUBLE_CELL(20000)
    measured=s.utils.support_types.SPICEDOUBLE_CELL(20000)
    files=[]
    for name in FILES:
        p=root/"Kernels/Cassini/ck"/name
        w=s.ckcov(str(p),-82000,False,"INTERVAL",0,"TDB")
        all_windows=s.wnunid(all_windows,w)
        if name.endswith("ra.bc"): measured=s.wnunid(measured,w)
        s.furnsh(str(p))
        files.append(dict(path=str(p.relative_to(root)), url="https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/ck/"+name,
                          sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size,
                          type="reconstructed" if name.endswith("ra.bc") else "predicted_from_as_flown_sequence",
                          intervals=s.wncard(w)))
    def intervals(w): return [s.wnfetd(w,i) for i in range(s.wncard(w))]
    def covered_seconds(w): return sum(max(0,min(end,b)-max(start,a)) for a,b in intervals(w))
    gaps=[]; cursor=start
    for a,b in intervals(all_windows):
        if b<start or a>end: continue
        if a>cursor: gaps.append((cursor,min(a,end)))
        cursor=max(cursor,b)
    if cursor<end: gaps.append((cursor,end))
    if any(b-a>10 for a,b in gaps): raise RuntimeError("Attitude gap exceeds 10 seconds")
    if not s.wnelmd(start,all_windows) or not s.wnelmd(end,all_windows):
        raise RuntimeError("Arc endpoints require unavailable attitude extrapolation")
    # Verify frame resolution at a one-minute grid and both gap boundaries.
    t=start
    while t<=end:
        if s.wnelmd(t,all_windows): s.pxform("CASSINI_SC_COORD","J2000",t)
        t+=60
    for a,b in gaps:
        s.pxform("CASSINI_SC_COORD","J2000",a);s.pxform("CASSINI_SC_COORD","J2000",b)
    result=dict(start_utc=args.start_utc,end_utc=s.et2utc(end,"ISOC",3),ck_id=-82000,
                frame="CASSINI_SC_COORD",reference="J2000",sclk="cas00172.tsc",files=files,
                reconstructed_seconds=covered_seconds(measured),
                as_flown_seconds=covered_seconds(all_windows)-covered_seconds(measured),
                gap_seconds=sum(b-a for a,b in gaps),max_bridge_seconds=10,
                gaps=[dict(start_utc=s.et2utc(a,"ISOC",6),end_utc=s.et2utc(b,"ISOC",6),seconds=b-a) for a,b in gaps],
                clock_note="Reconstructed CKs were produced with cas00086.tsc; encoded ticks are consistently mapped here with the later cas00172.tsc, also used to generate the as-flown CKs. UTC coverage is recomputed with this correlation.",
                policy="Reconstructed priority; predicted as-flown sequence only in telemetry gaps; SLERP only for bounded gaps <=10 s, no extrapolation.")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("files",)},indent=2))

if __name__=="__main__": main()
