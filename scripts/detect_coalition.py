#workspace
import graph_tool.all as gt
import scipy.stats as sc
import leidenalg as la
import igraph as ig
import pandas as pd
import numpy as np

#resolution
dl = []
op = la.Optimiser()
rs = np.round([0] + np.arange(0.85, 1.15, 0.01).tolist(), 2)
for r in rs:
    gr = ig.Graph.Read_GraphML("/m/triton/scratch/work/malkama5/congoCoalition/complete.xml")
    ws = np.array(gr.es["weight"])
    gr.es["leight"] = np.log(1 + ws)
    pts = []
    for seed in range(100):
        pt = la.RBConfigurationVertexPartition(gr, weights = gr.es["leight"], resolution_parameter = r)
        op.set_rng_seed(seed)
        op.optimise_partition(pt, n_iterations = -1)
        pts.append(pt.membership)
    gr = gt.load_graph("/m/triton/scratch/work/malkama5/congoCoalition/complete.xml")
    gr.ep.leight = gr.ep.weight.copy()
    gr.ep.leight.a = np.log(1 + gr.ep.leight.a)
    pm = gt.PartitionModeState(bs = pts, relabel = True, converge = True)
    pv = pm.get_marginal(gr)
    pr = []
    for v in gr.vertices():
        pr.append(np.argmax(pv[v]))
    b = gr.new_vp("int")
    for v in range(len(pr)):
        b[v] = pr[v]
    st = gt.BlockState(gr, b = b, recs = [gr.ep.leight], rec_types = ["real-exponential"])
    dl.append(st.entropy())
r = np.round(rs[dl.index(np.min(dl))], 2)
np.min(dl) - gt.BlockState(gr, B = 1, recs = [gr.ep.leight], rec_types = ["real-exponential"]).entropy()

#partition
gr = ig.Graph.Read_GraphML("/m/triton/scratch/work/malkama5/congoCoalition/complete.xml")
ws = np.array(gr.es["weight"])
gr.es["leight"] = np.log(1 + ws)
pts = []
for seed in range(1000):
    pt = la.RBConfigurationVertexPartition(gr, weights = gr.es["leight"], resolution_parameter = r)
    op.set_rng_seed(seed)
    op.optimise_partition(pt, n_iterations = -1)
    pts.append(pt.membership)
gr = gt.load_graph("/m/triton/scratch/work/malkama5/congoCoalition/complete.xml")
gr.ep.leight = gr.ep.weight.copy()
gr.ep.leight.a = np.log(1 + gr.ep.leight.a)
pm = gt.PartitionModeState(bs = pts, relabel = True, converge = True)
pv = pm.get_marginal(gr)
pr = []
for v in gr.vertices():
    pr.append(np.argmax(pv[v]))
b = gr.new_vp("int")
for v in range(len(pr)):
    b[v] = pr[v]
st = gt.BlockState(gr, b = b, recs = [gr.ep.leight], rec_types = ["real-exponential"])
gr.vp.block = b

#refine
sd = 1; gt.seed_rng(sd); np.random.seed(sd)
ps = gt.sfdp_layout(gr, groups = gr.vp.block, gamma = 2)
st.draw(pos = ps, vertex_shape = "pie", vertex_pie_fractions = pv, vertex_text = gr.vp.name, vertex_font_size = 6)

bx = len(np.unique(gr.vp.block.a))
for v in gr.vertices():
    if gr.vp.name[v] == "CD111":
        gr.vp.block[v] = bx
    if gr.vp.name[v] == "CD112":
        gr.vp.block[v] = bx
    if gr.vp.name[v] == "CD113":
        gr.vp.block[v] = bx

#coordinates
sd = 1; gt.seed_rng(sd); np.random.seed(sd)
ps = gt.sfdp_layout(gr, groups = gr.vp.block, gamma = 3)
#st.draw(pos = ps, vertex_color = "#000000", vertex_fill_color = gr.vp.block, vertex_text = gr.vp.block, vertex_font_size = 6)
for v in gr.vertices():
    ps[v][0] = -ps[v][0]
    ps[v][1] = -ps[v][1]
#st.draw(pos = ps, vertex_color = "#000000", vertex_fill_color = gr.vp.block, vertex_text = gr.vp.block, vertex_font_size = 6)
st.draw(pos = ps, vertex_color = "#000000", vertex_fill_color = gr.vp.block, vertex_text = gr.vp.name, vertex_font_size = 6)
for v in gr.vertices():
    if gr.vp.name[v] == "CD112":
        ps[v][0] = ps[v][0] - 7
        ps[v][1] = ps[v][1] - 8
    if gr.vp.name[v] == "CD113":
        ps[v][0] = ps[v][0] - 7
        ps[v][1] = ps[v][1] - 4
    if gr.vp.name[v] == "CD065":
        ps[v][0] = ps[v][0] + 9
        ps[v][1] = ps[v][1] + 4
    if gr.vp.name[v] == "CD172":
        ps[v][0] = ps[v][0] + 0
        ps[v][1] = ps[v][1] + 3
    if gr.vp.name[v] == "CD162":
        ps[v][0] = ps[v][0] + 1
        ps[v][1] = ps[v][1] - 2
    if gr.vp.name[v] == "CD143":
        ps[v][0] = ps[v][0] - 9
        ps[v][1] = ps[v][1] + 0
    if gr.vp.name[v] == "CD136":
        ps[v][0] = ps[v][0] + 9
        ps[v][1] = ps[v][1] + 1
    if gr.vp.name[v] == "CD053":
        ps[v][0] = ps[v][0] + 7
        ps[v][1] = ps[v][1] + 1
    if gr.vp.name[v] == "CD203":
        ps[v][0] = ps[v][0] + 2
        ps[v][1] = ps[v][1] - 4
st.draw(pos = ps, vertex_color = "#000000", vertex_fill_color = gr.vp.block, vertex_text = gr.vp.name, vertex_font_size = 6)

#attributes
at = pd.DataFrame({"index": gr.vertex_index, "actor": gr.vp.name, "block": gr.vp.block, "x": np.nan, "y": np.nan, "degree": np.nan, "wegree": np.nan, "closeness": gt.closeness(gr).a.tolist(), "betweenness": gt.betweenness(gr)[0].a.tolist()})

x = []
y = []
for v in gr.vertices():
    x.append(ps[v][0])
    y.append(ps[v][1])
at["x"] = x
at["y"] = y

nd = []
wd = []
for v in gr.vertices():
    nd.append(int(v.out_degree(weight = None)))
    wd.append(int(v.out_degree(weight = gr.ep.weight)))
at["degree"] = nd
at["wegree"] = wd

vc = at["block"].value_counts()
mp = {v:k for k, v in dict(enumerate(vc.index)).items()}
at["block"] = at["block"].map(mp)
at = pd.merge(at, pd.read_csv("/m/triton/scratch/work/malkama5/congoCoalition/attributes.csv", header = 0, sep = ";").loc[:, ["actor", "realm", "type", "name"]], on = "actor", how = "left")

at = at.sort_values(["block", "degree"])
at["coalition"] = np.nan
at.loc[at["block"] == 0, "coalition"] = "establishment+"
at.loc[at["block"] == 1, "coalition"] = "science_policy_regime"
at.loc[at["block"] == 2, "coalition"] = "forest_watchdog"
at.loc[at["block"] == 3, "coalition"] = "development_industry"
at.loc[at["block"] == 4, "coalition"] = "financial_justice"
at.loc[at["block"] == 5, "coalition"] = "family_farming"
at["order"] = np.arange(gr.num_vertices())
at = at.sort_values("index")
hdr = at.columns.tolist()
hdr = ";".join([str(e) for e in hdr])
np.savetxt("/m/triton/scratch/work/malkama5/congoCoalition/attributex.csv", at.values, delimiter = ";", fmt = "%s", header = hdr, comments = "")

#draw
cpie = ["#008E97FF", "#FC4C02FF", "#7FFFD4FF", "#151716FF", "#552583FF", "#FDB927FF", "#005778FF"]
vo = gr.new_vp("int")
cf = gr.new_vp("string")
ev = gr.new_vp("float")
ps = gr.new_vp("vector<double>")
for v in np.arange(gr.num_vertices()):
    ps[v] = [at["x"][v], at["y"][v]]
    ps[v][0] = ps[v][0] * 3
    ps[v][1] = ps[v][1] * 1
    ev[v] = at["degree"][v]
    if at["block"][v] == 0:
        cf[v] = cpie[0]
    if at["block"][v] == 1:
        cf[v] = cpie[1]
    if at["block"][v] == 2:
        cf[v] = cpie[2]
    if at["block"][v] == 3:
        cf[v] = cpie[3]
    if at["block"][v] == 4:
        cf[v] = cpie[4]
    if at["block"][v] == 5:
        cf[v] = cpie[5]
    if at["realm"][v] == "national":
        ps[v][0] = ps[v][0] + 10
        ps[v][1] = ps[v][1] + 60
        vo[v] = 0
    else:
        vo[v] = 1
eo = gr.new_ep("int")
for e in gr.edges():
    src, trg = e
    cmp = []
    for v in range(gr.num_vertices()):
        if gr.vertex_index[v] == src:
            cmp.append(at["realm"][v])
        if gr.vertex_index[v] == trg:
            cmp.append(at["realm"][v])
    if (cmp[0] == "national") & (cmp[1] == "national"):
        eo[e] = 0
    if (cmp[0] == "national") & (cmp[1] != "national"):
        eo[e] = 1
    if (cmp[0] != "national") & (cmp[1] == "national"):
        eo[e] = 1
    if (cmp[0] != "national") & (cmp[1] != "national"):
        eo[e] = 2
st.draw(pos = ps, vertex_fill_color = cf, vertex_color = cf, edge_pen_width = gt.prop_to_size(gr.ep.leight, 0.1, 0.3, log = False), vertex_size = gt.prop_to_size(ev, 2, 4, log = False), vertex_aspect = 3, vorder = vo, eorder = eo, output = "/m/triton/scratch/work/malkama5/congoCoalition/blocks.svg")

st.draw(pos = ps, vertex_shape = "pie", vertex_pie_fractions = pv, edge_pen_width = gt.prop_to_size(gr.ep.leight, 0.1, 0.3, log = False), vertex_size = gt.prop_to_size(ev, 4, 8, log = False), vorder = vo, eorder = eo, output = "/m/triton/scratch/work/malkama5/congoCoalition/blocksFractions.svg")

#nested
sn = gt.NestedBlockState(gr, bs = [gr.vp.block, gr.vp.block])
od = gr.new_vp("int")
for v in np.arange(gr.num_vertices()):
    od[v] = at["order"][v]
sd = 1; gt.seed_rng(sd); np.random.seed(sd)
gt.draw_hierarchy(sn, rel_order = od, vertex_fill_color = cf, vertex_color = cf, eorder = gr.ep.leight, edge_pen_width = gt.prop_to_size(gr.ep.leight, 0.2, 0.6, log = False), vertex_size = 8, hvertex_size = 0, hedge_pen_width = 0, hedge_marker_size = 0, output = "/m/triton/scratch/work/malkama5/congoCoalition/nestedTopicAll.svg")[0]

cfi = []
for v in gr.vertices():
    cfi.append(cf[v])
at["cx"] = cfi

#topics
gs = [
    gt.load_graph("/m/triton/scratch/work/malkama5/congoCoalition/complete.xml"),
    gt.load_graph("/m/triton/scratch/work/malkama5/congoCoalition/agriculture.xml"),
    gt.load_graph("/m/triton/scratch/work/malkama5/congoCoalition/bioenergy.xml"),
    gt.load_graph("/m/triton/scratch/work/malkama5/congoCoalition/economy.xml"),
    gt.load_graph("/m/triton/scratch/work/malkama5/congoCoalition/environment.xml"),
    gt.load_graph("/m/triton/scratch/work/malkama5/congoCoalition/equity.xml"),
    gt.load_graph("/m/triton/scratch/work/malkama5/congoCoalition/extractive.xml"),
    gt.load_graph("/m/triton/scratch/work/malkama5/congoCoalition/forest.xml"),
    gt.load_graph("/m/triton/scratch/work/malkama5/congoCoalition/policy.xml")]
for i, g in enumerate(gs):
    ga = []
    for v in g.vertices():
        ga.append(g.vp.name[v])
    a = np.isin(at["actor"], ga)
    a = at["actor"][a == False].tolist()
    g.add_vertex(len(a))
    a = ga + a
    for index, vertex in enumerate(g.vertices()):
        g.vp.name[vertex] = a[index]    
    a = pd.merge(pd.DataFrame({"actor": g.vp.name}), at, on = "actor", how = "left")
    cxg = g.new_vp("string")
    blg = g.new_vp("int")
    odg = g.new_vp("int")
    g.ep.leight = g.ep.weight.copy()
    g.ep.leight.a = np.log(1 + g.ep.leight.a)
    for v in range(g.num_vertices()):
        cxg[v] = a["cx"][v]
        blg[v] = a["block"][v]
        odg[v] = a["order"][v]
    for v in g.vertices():
        if v.out_degree() == 0:
            cxg[v] = "#FFFFFF00"
    n = gt.NestedBlockState(g, bs = [blg, blg])
    gt.draw_hierarchy(n, rel_order = odg, vertex_fill_color = cxg, vertex_color = cxg, eorder = g.ep.leight, edge_pen_width = gt.prop_to_size(g.ep.leight, 0.2, 0.6, log = False), vertex_size = 8, hvertex_size = 0, hedge_pen_width = 0, hedge_marker_size = 0, output = "/m/triton/scratch/work/malkama5/congoCoalition/nestedTopic" + str(i) + ".svg")

#exin
def exin(elist, vlist, weights = False, adaptive = False, entire = False):
    """External-Internal Index by Krackhardt and Stern (1988) based on pairwise comparisons between blocks. Takes two Pandas dataframes as input. The 'elist' needs to contain all edges with 'source' and 'target' columns; if 'weights' is True, an additional column 'weight' must be present. The 'vlist' needs to contain all vertices with 'vertex' and 'member' columns. If 'adaptive' is True, the index is normalised by the product of the block sizes. If 'entire' is False, a block matrix with block sizes along the diagonal is returned. If 'entire' is True, the weighted average of the pairwise comparisons is returned with weights corresponding to the product of the block sizes. Requires Numpy 'as np' and Pandas 'as pd'."""

    el = elist
    nl = vlist
    nl = nl.rename(columns = {"vertex": "source", "member": "member_s"})
    el = pd.merge(el, nl, on = "source", how = "left")
    nl = nl.rename(columns = {"source": "target", "member_s": "member_t"})
    el = pd.merge(el, nl, on = "target", how = "left")
    nl = nl.rename(columns = {"target": "vertex", "member_t": "member"})
    el["internal"] = np.where(el["member_s"] == el["member_t"], 1, 0)
    el["external"] = np.where(el["member_s"] != el["member_t"], 1, 0)
    
    if weights:
        el["internal_w"] = np.where(el["internal"] == 1, el["weight"], 0)
        el["external_w"] = np.where(el["external"] == 1, el["weight"], 0)
        el["internal"] = el["internal_w"]
        el["external"] = el["external_w"]
    
    bname, bsize = np.unique(nl["member"], return_counts = True)
    bs = pd.DataFrame(np.zeros(shape = (len(bsize), len(bsize))))
    for j in np.flip(np.arange(len(bname))):
        bs = bs.rename(columns = {j: bname[j]}, index = {j: bname[j]})
    np.fill_diagonal(bs.values, bsize)
    
    bm = bs.copy()
    np.fill_diagonal(bm.values, 0)
    for edge in range(len(el)):
        s = el["member_s"][edge]
        t = el["member_t"][edge]
        if s == t:
            bm.loc[s, t] += el["internal"][edge]
        else:
            bm.loc[s, t] += el["external"][edge]
    bm = bm.add(bm.transpose())
    np.fill_diagonal(bm.values, np.diag(bm) / 2)

    a = bname.tolist() * len(bname)
    b = np.repeat(bname, len(bname)).tolist()

    if adaptive:
        bn_internal = bs.copy()
        np.fill_diagonal(bn_internal.values, 0)
        for j in bname:
            bn_internal.loc[j, j] = (bm.loc[j, j]) / (bs.loc[j, j] * (bs.loc[j, j] - 1) * 0.5)

        bs_divider = bs.copy()
        np.fill_diagonal(bs_divider.values, 0)
        for k in a:
            for l in b:
                bs_divider.loc[k, l] = (bs.loc[k, k] * (bs.loc[l, l] - 1) * 0.5)

        bn = bm.div(bs_divider)
        np.fill_diagonal(bn.values, np.diag(bn_internal))
    else:
        bn = bm.copy()
    
    ei = bs.copy()
    np.fill_diagonal(ei.values, 0)
    for k in a:
        for l in b:
            lefthand = -(bn.loc[k, k] + bn.loc[l, l] - bn.loc[k, l] - bn.loc[l, k])
            righthand = (bn.loc[k, k] + bn.loc[l, l] + bn.loc[k, l] + bn.loc[l, k])
            if righthand != 0:
                ei.loc[k, l] = lefthand / righthand
            else:
                ei.loc[k, l] = 0
    np.fill_diagonal(ei.values, np.diag(bs))

    if entire:
        bs_product = bs.copy()
        np.fill_diagonal(bs_product.values, 0)
        for k in a:
            for l in b:
                bs_product.loc[k, l] = bs.loc[k, k] * bs.loc[l, l]
        a_values = ei.values[np.triu_indices_from(ei.values, k = 1)]
        a_weights = bs_product.values[np.triu_indices_from(bs_product.values, k = 1)]
        if sum(a_weights) == 0:
            ei = 0
        else:
            ei = np.average(a_values, weights = a_weights)

    return ei

#polarisation
ex = []
for g in gs:
    src = []
    trg = []
    wgh = []
    for e in g.edges():
        sc, tg = e
        src.append(str(sc))
        trg.append(str(tg))
        wgh.append(g.ep.weight[e])
    el = pd.DataFrame({"source": src, "target": trg, "weight": wgh})
    vd = []
    for v in g.vertices():
        vd.append(str(g.vertex_index[v]))
    vl = pd.merge(pd.DataFrame({"vertex": vd, "actor": g.vp.name}), at.loc[:, ["actor", "block"]], on = "actor", how = "left")
    vl = vl.rename(columns = {"block": "member"})
    ex.append(exin(elist = el, vlist = vl, weights = True, adaptive = True, entire = False))
for e in ex:
    e.index = ["establishment+", "science_policy_regime", "forest_watchdog", "development_industry", "financial_justice", "family_farming"]
    e.columns = ["establishment+", "science_policy_regime", "forest_watchdog", "development_industry", "financial_justice", "family_farming"]
np.round(ex[0] * - 1, 2)

#mentions
gm = gt.load_graph("/m/triton/scratch/work/malkama5/congoCoalition/mn.xml")
ga = pd.merge(pd.DataFrame({"actor": gm.vp.name}), pd.read_csv("/m/triton/scratch/work/malkama5/congoCoalition/attributes.csv", header = 0, sep = ";").loc[:, ["actor", "realm"]], on = "actor", how = "left")
ly = gm.new_vp("int")
for v in range(gm.num_vertices()):
    if ga["realm"][v] != "national":
        ly[v] = 1
    if ga["realm"][v] == "national":
        ly[v] = 2
sd = 6; gt.seed_rng(sd); np.random.seed(sd)
ps = gt.sfdp_layout(gm, groups = ly, gamma = 3)
st = gt.BlockState(gm, b = ly)
for v in gm.vertices():
    if v.out_degree() == 0:
        ps[v][0] = ps[v][0] * 0.3
        ps[v][1] = ps[v][1] * 0.3
#st.draw(pos = ps, vertex_text = ly, edge_end_marker = "arrow", edge_marker_size = 10, edge_control_points = [])
#st.draw(pos = ps, vertex_text = gm.vp.name, edge_end_marker = "arrow", edge_marker_size = 10, edge_control_points = [])
for v in gm.vertices():
    if gm.vp.name[v] == "CD178":
        ps[v][0] = ps[v][0] - 6
        ps[v][1] = ps[v][1] + 0
    if gm.vp.name[v] == "CD022":
        ps[v][0] = ps[v][0] - 7
        ps[v][1] = ps[v][1] - 7
    if gm.vp.name[v] == "CD201":
        ps[v][0] = ps[v][0] + 1
        ps[v][1] = ps[v][1] - 2
    if gm.vp.name[v] == "CD052":
        ps[v][0] = ps[v][0] + 1
        ps[v][1] = ps[v][1] - 2
    if gm.vp.name[v] == "CD111":
        ps[v][0] = ps[v][0] - 2
        ps[v][1] = ps[v][1] - 5
    if gm.vp.name[v] == "CD112":
        ps[v][0] = ps[v][0] - 3
        ps[v][1] = ps[v][1] - 14
    if gm.vp.name[v] == "CD203":
        ps[v][0] = ps[v][0] - 0
        ps[v][1] = ps[v][1] - 4
    if gm.vp.name[v] == "CD193":
        ps[v][0] = ps[v][0] - 6
        ps[v][1] = ps[v][1] + 3
    if gm.vp.name[v] == "CD189":
        ps[v][0] = ps[v][0] - 5
        ps[v][1] = ps[v][1] + 2
    if gm.vp.name[v] == "CD059":
        ps[v][0] = ps[v][0] - 4
        ps[v][1] = ps[v][1] + 0
    if gm.vp.name[v] == "CD110":
        ps[v][0] = ps[v][0] - 4
        ps[v][1] = ps[v][1] + 0
    if gm.vp.name[v] == "CD004":
        ps[v][0] = ps[v][0] - 3
        ps[v][1] = ps[v][1] + 0
    if gm.vp.name[v] == "CD158":
        ps[v][0] = ps[v][0] + 4
        ps[v][1] = ps[v][1] + 1
    if gm.vp.name[v] == "CD194":
        ps[v][0] = ps[v][0] + 1
        ps[v][1] = ps[v][1] + 4
    if gm.vp.name[v] == "CD048":
        ps[v][0] = ps[v][0] + 2
        ps[v][1] = ps[v][1] + 0
    if gm.vp.name[v] == "CD164":
        ps[v][0] = ps[v][0] + 2
        ps[v][1] = ps[v][1] + 0
    if gm.vp.name[v] == "CD068":
        ps[v][0] = ps[v][0] + 1
        ps[v][1] = ps[v][1] + 0
    if gm.vp.name[v] == "CD080":
        ps[v][0] = ps[v][0] - 2
        ps[v][1] = ps[v][1] + 2
    if gm.vp.name[v] == "CD053":
        ps[v][0] = ps[v][0] + 0
        ps[v][1] = ps[v][1] + 2
    if gm.vp.name[v] == "CD142":
        ps[v][0] = ps[v][0] + 0
        ps[v][1] = ps[v][1] + 3
    if gm.vp.name[v] == "CD131":
        ps[v][0] = ps[v][0] + 0
        ps[v][1] = ps[v][1] - 3
    if gm.vp.name[v] == "CD124":
        ps[v][0] = ps[v][0] + 3
        ps[v][1] = ps[v][1] + 0
    if gm.vp.name[v] == "CD122":
        ps[v][0] = ps[v][0] + 3
        ps[v][1] = ps[v][1] + 0
    if gm.vp.name[v] == "CD143":
        ps[v][0] = ps[v][0] + 2
        ps[v][1] = ps[v][1] + 0
    if gm.vp.name[v] == "CD149":
        ps[v][0] = ps[v][0] + 2
        ps[v][1] = ps[v][1] + 0
    if gm.vp.name[v] == "CD139":
        ps[v][0] = ps[v][0] + 0
        ps[v][1] = ps[v][1] + 2
    if gm.vp.name[v] == "CD160":
        ps[v][0] = ps[v][0] + 0
        ps[v][1] = ps[v][1] + 2
    if gm.vp.name[v] == "CD151":
        ps[v][0] = ps[v][0] + 2
        ps[v][1] = ps[v][1] + 0
st.draw(pos = ps, vertex_text = ly, edge_end_marker = "arrow", edge_marker_size = 10, edge_control_points = [])
au = gt.hits(gm)[1]
cf = gm.new_vp("string")
vo = gm.new_vp("int")
for v in gm.vertices():
    if st.b == np.unique(st.b.a)[0]:
        ps[v][0] = ps[v][0] + 0
        ps[v][1] = ps[v][1] + 2
    if st.b == np.unique(st.b.a)[1]:
        ps[v][0] = ps[v][0] + 0
        ps[v][1] = ps[v][1] - 2
    if gm.vp.layership[v] == 1:
        ps[v][0] = ps[v][0] - 5
        ps[v][1] = ps[v][1] - 30
        vo[v] = 0
        cf[v] = "#007fffcc"
    if gm.vp.layership[v] == 2:
        ps[v][0] = ps[v][0] + 0
        ps[v][1] = ps[v][1] + 0
        vo[v] = 1
        cf[v] = "#ce102fcc"
    ps[v][0] = ps[v][0] * 3
    ps[v][1] = ps[v][1] * 1
eo = gm.new_ep("int")
for e in gm.edges():
    src, trg = e
    cmp = []
    for v in range(gm.num_vertices()):
        if gm.vertex_index[v] == src:
            cmp.append(ga["realm"][v])
        if gm.vertex_index[v] == trg:
            cmp.append(ga["realm"][v])
    if (cmp[0] == "national") & (cmp[1] == "national"):
        eo[e] = 0
    if (cmp[0] == "national") & (cmp[1] != "national"):
        eo[e] = 1
    if (cmp[0] != "national") & (cmp[1] == "national"):
        eo[e] = 1
    if (cmp[0] != "national") & (cmp[1] != "national"):
        eo[e] = 2
st.draw(pos = ps, vertex_fill_color = cf, vertex_color = "#f7d6184d", vertex_size = gt.prop_to_size(au, 2, 4, log = False), vertex_aspect = 3, vorder = vo, eorder = eo, edge_pen_width = 0.1, edge_end_marker = "arrow", edge_marker_size = 0, edge_control_points = [], output = "/m/triton/scratch/work/malkama5/congoCoalition/mentions.svg")

