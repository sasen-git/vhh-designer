import os, glob, pickle

def find_pkls(root):
    return sorted(glob.glob(os.path.join(root, "**", "analysis_epistasis_legacy.pkl"), recursive=True))

def summarize(pkl_path):
    ep = pickle.load(open(pkl_path, "rb"))
    clusters = ep.get("analysis_2_vernier_clusters", {})
    total = len(clusters)

    singletons = n10 = n100 = 0
    clustered_n = 0

    for v in clusters.values():
        n = int(v.get("cdr3_length", {}).get("n", 0))
        clustered_n += n
        if n == 1: singletons += 1
        if n >= 10: n10 += 1
        if n >= 100: n100 += 1

    singleton_pct = 100 * singletons / total if total else 0
    return total, singleton_pct, n10, n100, clustered_n


def run(root):
    print(f"\n=== SUMMARY: {root} ===")
    print("run\tclusters\tsingleton%\tn>=10\tn>=100\tclustered_n")
    for pkl in find_pkls(root):
        name = os.path.basename(os.path.dirname(os.path.dirname(pkl)))
        total, sp, n10, n100, cn = summarize(pkl)
        print(f"{name}\t{total:,}\t{sp:.1f}\t{n10:,}\t{n100:,}\t{cn:,}")


if __name__ == "__main__":
    import sys
    for root in sys.argv[1:]:
        run(root)
