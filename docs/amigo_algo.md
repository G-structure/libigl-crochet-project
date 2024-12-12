# AmiGo: Computational Design of Amigurumi Crochet Patterns

## 2. REPRESENTATION

### 2.1 Background and Notations
Given a closed manifold triangle mesh \( M = (V, E) \), a seed vertex \( s \in V \), and a stitch width \( w \in \mathbb{R} \), our goal is to generate human-readable instructions \( P(M, s, w) \) for crocheting \( M \) from the point \( s \), with the given stitch width \( w \).

Crochet has a wide variety of stitches, and we focus here on the simple stitch used for Amigurumi, named single crochet (sc). This is an approximately square stitch, thus covering \( M \) with sc stitches is equivalent to constructing a quad re-mesh of \( M \), where each quad is a square, and all the edge lengths are constant. This is of course not possible unless the surface is developable, i.e., has zero Gaussian curvature. In practice, curved geometry is accommodated in crochet by introducing stitches which locally increase \( inc(x) \) or decrease \( dec(x) \) the amount of stitches by \( x \). Figure 3 (left) shows an example of the \( inc \) and \( dec \) stitches on a crocheted patch. Crochet instructions for Amigurumi typically include rows, where each row is a series of \( sc \), \( inc \), \( dec \) stitches. Figure 2 (middle) shows the instructions (pattern) for crocheting the sphere in Figure 2 (right).

### 2.2 The Crochet Graph
A crochet stitch is composed of a top, a base, and a stem, where the \( inc \), \( dec \) stitches have multiple stems, see Figure 3 (right). The top of one stitch is always the base of some stitch on the next row, and similarly, each stitch has a base on the previous row. Therefore, a natural abstraction of the stitch pattern is to consider the stitches and their interconnections as a graph.

Specifically, we define the Crochet Graph \( G = (S, R \cup C) \), whose vertices \( S \) are tops/bases of stitches, where a vertex \( (i, j) \in S \) is the base of the \( j \)-th stitch in row \( i \), and the vertices in each row are consecutively ordered. The column edges \( C \) are stems of stitches, and the connectivity between the bases in each row is represented by the row edges \( R \). We denote the total number of rows by \( N \). Figure 2 (left) shows the crochet graph corresponding to the crocheted sphere in Figure 2 (right).

A crochet graph is an intermediate representation between the input triangle mesh \( M \), and the output instructions \( P \). Our goal is to generate a graph such that it is (1) translatable to valid crochet instructions \( P \), and (2) when \( P \) is crocheted and stuffed, the result resembles the input mesh. Note that there exist multiple instructions \( P \) for the same graph \( G \), and within this space we aim for instructions which are human-readable.

We base our algorithm on the following observations:

**Definition 2.1.** A coupling [Gold and Sharir 2018] \( C = (c_1, \ldots, c_k) \) between two sequences \( A = (p_1, \ldots, p_n) \) and \( B = (q_1, \ldots, q_m) \) is an ordered sequence of distinct pairs of points from \( A \times B \), such that \( c_1 = (p_1, q_1) \), \( c_k = (p_n, q_m) \) and 
\[
c_r = (p_s, q_t) \implies c_{r+1} \in \{ (p_{s+1}, q_t), (p_s, q_{t+1}), (p_{s+1}, q_{t+1}) \}, \forall r < k.
\]

**Definition 2.2.** Let \( S_i, S_{i+1}, 1 \leq i < N \), be the vertices of two consecutive rows of \( G = (S, R \cup C) \), where \( S_i = \{ (i, 1), \ldots, (i, n_i) \} \), where \( (i, j) \in S \), and \( n_i \) is the number of vertices in row \( i \). If there exists a coupling \( C \) between \( S_i \) and \( S_{i+1} \) such that for all \( p_s \in S_i, q_t \in S_{i+1} \), we have that \( (p_s, q_t) \in C \) if and only if \( (p_s, q_t) \in C \), then the two rows are coupled.

**Observation 2.3.** If all the pairs of consecutive rows of \( G \) are coupled, then there exist valid crochet instructions \( P(G) \) that use only the instructions \( sc \), \( inc(x) \), and \( dec(x) \).

**Definition 2.4.** Let \( X_G : S \to M \) be an embedding of the vertices of \( G \) on \( M \). An embedded edge of \( G \) is a shortest geodesic between the embedding of two vertices of \( G \) which share an edge, or between the embedding of the first and last vertices on the same row.

**Definition 2.5.** Let \( (p, q) \in M \) be two points whose geodesic distance is larger than some constant that depends on the stitch width \( w \), and let \( \gamma_{p,q} \) be the shortest geodesic between them. If \( \gamma_{p,q} \) intersects some embedded edge of an embedding \( X_G \), for any two such points, then we say that \( X_G \) covers \( M \).

**Observation 2.6.** Let \( X_G : S \to M \) be an embedding of the vertices of \( G \) which covers \( M \), and \( P(G) \) valid crochet instructions for \( G \). If all the edge lengths induced by \( X_G \) are equal to \( w \), then when \( P(G) \) is crocheted and stuffed the result will be "similar" to \( M \).

We discuss **Observation 2.3** in Section 5.1, where we show how to translate the graph into valid instructions. The **Observation 2.6** is in fact true only for a subset of meshes, as we discuss in the next section.

### 2.3 Crochetable Models
**Curvature.** Crocheting the patterns yields an empty flexible shell of fabric, which obtains its final shape by stuffing. Whether the stuffed model obtains the intended shape depends on how the model is stuffed (lightly, firmly, uniformly), as the yarn has some flexibility and will extend to accommodate if the model is over-stuffed. We will assume that the model is stuffed enough to attain maximal volume, but not too much to cause the yarn to stretch and generate gaps. Thus, we expect the resulting stitch size to be similar to the edge lengths induced by the embedding of the crochet graph.

If for a given graph \( G \), its embedding in 3D with edge lengths \( w \) is unique, we expect the crocheted and stuffed shape to be similar to the input surface.

Importantly, unless the shape is convex, the edge lengths alone (i.e., the metric) do not contain enough information to uniquely determine the shape of a non-stuffed model. For example, a surface that has a "crater" leads to edge lengths which can be realized either as a crater or as a "hill." However, if we add the maximal volume assumption, only the hill is possible. This in fact implies that surfaces which have "craters," or more formally, regions of negative mean curvature with positive Gaussian curvature, cannot be realized by crocheting and stuffing alone.

This is similar to the observation made by Konaković et al. [2018], that surfaces which have negative mean curvature cannot be realized by maximizing the volume with a given conformal metric (i.e., only isotropic scaling is allowed relative to the flat configuration). We handle this case similarly, by preprocessing the model so that it does not contain "craters" (Section 7.1.1).

Furthermore, in our case, since we allow anisotropic scaling, negative mean curvature with negative Gaussian curvature is in fact possible, but requires a modified sampling rate, which we discuss in Section 7.1.2. To conclude, in terms of geometric obstructions to crochetability, the four possible curvature situations are summarized in Table 1. Note that, as with any sampling-dependent method, if the sampling rate (namely, the number of rows \(N\)) is too small compared to the feature size of the model, the crocheted output will lose some of the geometric detail.

**Branching.** Observation 2.6 requires that the crochet graph embedding \(X_G\) covers the input surface \(M\). Because of the special structure of this graph, this induces additional constraints on the possible geometries. Intuitively, models that branch (see Figure 7) cannot be covered in this way. Mathematically, this means that the geodesic distance function on \(M\) from the embedding of the seed vertex \(s\) cannot have saddles. This is solved by segmenting the shape, and crocheting the segments in an iterative manner. We explain this in detail in Section 7.2.

We first explain the generation of the crochet pattern for a simple non-branching model with positive mean curvature, and then discuss how we handle negative mean curvature and branching.

## 3. OVERVIEW
Given a 3D mesh \(M\), a seed point \(s\), and a stitch width \(w\), we first compute a crochet graph \(G\) and its embedding \(X_G\) such that they adhere to Observations 2.3 and 2.6 (Section 4). Then we compute the crochet pattern \(P(G)\) (Section 5).

To generate \(G\) and \(X_G\), we first compute the embedding of the vertices \(S\) on \(M\) (Section 4.1), and then derive from that the connectivity of \(G\), i.e., the row edges \(R\) and column edges \(C\) (Section 4.2).

To compute the pattern \(P(G)\), we first translate the graph into a program using standard code synthesis tools (Section 5.1), and then apply loop unrolling to make the pattern human-readable (Section 5.2). See Algorithm 1.

---

**Algorithm 1: An outline of our algorithm**  
**Input:** A triangle mesh \(M\), seed \(s\), stitch width \(w\)  
**Output:** Embedded crochet graph \(G = (S, R \cup C)\), \(X_G\), crochet pattern \(P(G)\)  

1. **Mesh to Graph** // Section 4  
    - **Geometry \(S, X_G\)** // Section 4.1  
    - **Connectivity \(R, C\)** // Section 4.2  
2. **Graph to Pattern** // Section 5  
    - **Graph to Program** // Section 5.1  
    - **Program to Pattern \(P(G)\)** // Section 5.2  

---

## 4. MESH TO CROCHET GRAPH

### 4.1 Geometry
Observation 2.3 implies that the vertices \(S\) should be grouped into ordered rows, where in each row the vertices have a well-defined order. We address this requirement by computing two non-negative monotonically increasing, constant speed functions \(f, g : M \to \mathbb{R}\) which define the row-order and column-order of every point on \(M\). Furthermore, Observation 2.6 implies that the distance between embedded rows, and between embedded vertices in the same row should be \(w\). We address this by sampling \(f, g\) appropriately.

**Row order \(f\).** Our models are closed (so they can be stuffed), and therefore the first and last rows in the graph \(G\) contain a single vertex. The first row contains only the seed \(s\), and its function value is \(f(s) = 0\). We take \(f(v), v \in V\) to be \(f(v) = d(v, s)\), where \(d\) is the geodesic distance. Thus, the isolines of \(f\) are rows, and two points \(p, q \in M\) are on the same row if \(f(p) = f(q)\). If \(f\) has more than one maximum, then we need to handle branching (see Section 7.2). Otherwise, the vertex that attains the maximum of \(f\), denoted as \(f_M\), will be the single vertex on the last row.

**Column order \(g\).** We first cut \(M\) along a geodesic from \(s\) to \(f_M\), so that our model and the graph that we compute have the same topology, and denote the cut model by \(M_C\). The requirements are that within each row the vertices of \(G\) have a well-defined order. A row is an isoline of \(f\), and therefore the rate of change along the isoline is given by the directional derivative of \(g\) in the direction of the tangent to the isoline. Specifically, the tangent to the isoline of \(f\) at a point \(p \in M\) is given by \(J \nabla f\), where \(J\) is the rotation by \(\pi/2\) in the tangent plane of \(p\). Thus, to find \(g\), we solve an optimization problem whose objective is to minimize:
\[
\int_{M_C} | \langle J \nabla f, \nabla g \rangle - 1 |^2, \quad \text{s.t. } g(B) = 0.
\]
Here, \(B \subset V\) is the longest connected path of boundary vertices of \(M_C\) along which \(f\) is strictly monotone.

**Sampling.** The functions \(f, g\) define a parameterization of \(M_C\) to the plane. We conjecture that this parameterization is bijective (as it was in all of our experiments), but leave the proof to future work. The parameterization may have a large metric distortion. However, if \(f(p) = f(q) = f_0\) for some two points \(p, q \in M\), then \(|g(p) - g(q)|\) is equal to the length of the isoline of \(f_0\) between \(p\) and \(q\). Therefore, we uniformly sample \(f, g\) on a \(2D\) grid of width \(w\), yielding the vertices of \(S\) with indices \((f/w, g/w)\). Pushing forward the sampled points to the mesh \(M_C\) yields the embedding of \(S\) on \(M_C\) (and therefore \(M\)), namely \(X_G\).

### 4.2 Connectivity

**Row edges \(R\):** Each two consecutive vertices of \(S\) on the same row are connected by a row edge. Namely, \(R = \bigcup_{i=1}^N R_i\), and \(R_i = \{((i, j), (i, j+1)) \,|\, j \in \{1, \ldots, n_i - 1\}\}\). Here \(n_i = |S_i| = |\{(i, j) \in S\}|\), namely the number of vertices in the \(i\)-th row.

Let \(x, y \in S\) be two consecutive vertices on the \(i\)-th row. Then we have that \(f(x) = f(y) = f_0\) and \(|g(x) - g(y)| = w\). Therefore, 
\[
d_\gamma(f_0)(X_G(x), X_G(y)) = w,
\]
where \(\gamma(f_0)\) is the isoline of \(f_0\) on \(M\), and \(d_\gamma(f_0)\) is the distance along the isoline. Hence, the Euclidean distance between the embedded vertices \(|X_G(x) - X_G(y)| \leq w\), and the distance tends to \(w\) for a "small enough" stitch size. Here, "small enough" means on the order of the square root of the radius of curvature of \(\gamma(f_0)\), which is given in terms of the normal curvature in direction \(J \nabla f\).

**Column edges \(C\):** First, Observation 2.3 requires that all pairs of consecutive rows are coupled. Let \(C_i\) be the coupling corresponding to rows \(S_i, S_{i+1}\), and let \((p_s, q_t) \in C_i\). Since \(p_s\) and \(q_t\) are on consecutive rows, and therefore embedded on isolines of \(f\) which differ by \(w\), the minimal distance \(|X_G(p_s) - X_G(q_t)|\) is close to \(w\). Therefore, if among all couplings we seek the minimizer of:
\[
\min_{C_i:\, \text{coupling}} \sum_{(p_s, q_t) \in C_i} |X_G(p_s) - X_G(q_t)|,
\]
then the length of the column edges will be close to \(w\).

A minimal coupling between every pair of consecutive rows is found by **Dynamic Time Warping (DTW)** [Gold and Sharir 2018; Sakoe and Chiba 1978].

---

## 5. CROCHET GRAPH TO INSTRUCTIONS

### 5.1 Graph to Program

In order to turn the crochet graph into instructions, we rely on the following observation: crochet instructions (patterns) constitute an instruction set, and as such, crocheting is an execution, and the finished object is an execution result. Moreover, because of the nature of a crocheted object, it is not only a result but a step-by-step execution trace.

Therefore, given a crochet object, or in this case its graph representation \(G\), deriving instructions constitutes a form of **execution reconstruction** [Zuo et al. 2021], a method for reconstructing the set of instructions that lead to an execution trace. While reconstructing an execution trace generally requires searching an exponential space of possible instruction sequences, the crochet instruction set is limited enough that reconstituting a trace is done in linear time using a **transducer**.

In order to reconstruct the trace, the degree of the graph vertices in both the current row \(S_i\) and next row \(S_{i+1}\) must be considered. Therefore, an execution reconstruction transducer accepts two consecutive rows and their connecting edges, and advances on the rows based on the vertex degrees in either row. For example:
- If at \(S_{i+1}\), the transducer is at a vertex with one edge connecting it to the previous row \(S_i\), and at the current head of \(S_i\) there is only a connection to the current head of \(S_{i+1}\), the transducer will yield a \(sc\).
- If there are \(x > 1\) connected vertices in \(S_i\) to the head of \(S_{i+1}\), the transducer will advance on \(S_i\) to consume them all, advance on a single vertex on \(S_{i+1}\), and produce a \(dec(x)\). 
- Analogously, for \(x > 1\) connected vertices in \(S_{i+1}\) to the head of \(S_i\), the transducer will advance on \(S_{i+1}\) to consume them, advance on a single vertex on \(S_i\), and then produce an \(inc(x)\). 

Figure 5 (left, middle) shows an illustration of the different situations. This is entirely analogous to the way crocheting the row is done.

Because the transducer chooses its transition based on which vertex at the current heads of the rows has a degree larger than 1, it is not technically deterministic. However, for a non-deterministic choice to exist, both the head of \(S_i\) should be connected to multiple vertices in \(S_{i+1}\) and vice versa. This is impossible since all the pairs of consecutive rows of \(G\) are coupled (see Figure 5 (right)). Thus, under the constraints of the input, i.e., a valid crochet graph \(G\), all the transitions are mutually exclusive, rendering the transducer’s behavior essentially deterministic.

### 5.2 Program to Human-readable Pattern

The instructions in a reconstituted trace can get quite repetitive. The raw output can contain the following rows:

row 2: sc, inc, sc, sc, sc, inc, sc, sc, sc, inc, sc, sc
row 3: sc, inc, sc, sc, sc, inc, sc, sc, sc, inc, sc, sc

indicating that row 2 is constructed via an \(sc\) instruction followed by an \(inc\), then another two \(sc\), repeating three times, then row 3 is constructed the same way. To make the instructions both succinct and more human-readable, it is customary to convert this code to:

rows 2-3: (sc, inc, 2sc)*3

In order to do this, we employ a disassembly technique called **loop folding** [Lee et al. 1994], which finds maximal repetitions of instructions and turns them into loops. Loop folding of crochet instructions occurs at three different levels: repeating rows, repeating sequences of stitches, and repeating stitches. 

In order to find maximal repeating sequences, the order in which loops are folded must be:
1. Sequences first,
2. Then repeating stitches.

For instance, in the above example, if repeating stitches were folded first, row 2 after this initial folding would be \(sc, inc, 3sc, inc, 3sc, inc, 2sc\), which means the repetition that can be identified within the row would have been smaller. Instead, it is first folded to find the repeating sequence \((sc, inc, sc, sc)*3\), which is maximal, and then the internal sequence is re-folded to identify the \(2sc\). Finally, identical rows are folded together.

## 6. CROCHET GRAPH TO 3D EMBEDDING

We generate a 3D embedding \(Y_G\) of the crochet graph \(G\) using **ShapeUp** [Bouaziz et al. 2012], in order to obtain a visualization of the expected result. Interestingly, using purely geometric conditions, the expected result is quite similar to the crocheted model in practice.

We use the following constraints for ShapeUp:
1. The column edges which represent \(sc\) stitches, as well as the row edges, are constrained to have length \(w\).
2. The embedding of the seed point \(s\) is fixed to the position of the seed vertex.
3. Smoothness.

We initialize \(Y_G\) using the sampled points \(X_G\). Some of the figures (see for example Figure 6) show, in addition to the embedded crochet graph \(X_G\), the ShapeUp result \(Y_G\).

---

## 7. OBSTRUCTIONS TO CROCHETABILITY

### 7.1 Negative Mean Curvature

As curvature computations are not scale-invariant, all the models are normalized to have a surface area equal to 1, so that we can use the same parameters across all the models.

#### 7.1.1 Positive Gaussian Curvature
Crocheting and stuffing a model which has "craters," i.e., regions of positive Gaussian curvature and negative mean curvature, will not yield a geometry that is similar to the input mesh. Thus, we apply a pre-processing step (similarly to Konaković et al. [2018]) to smooth out the craters. Specifically, we apply **Conformal Mean Curvature Flow** [Kazhdan et al. 2012] localized to these areas until the mean curvature is positive everywhere.

#### 7.1.2 Negative Gaussian Curvature
The sampling rate of the isolines of \(f\) is determined by the directional derivative of \(g\) with respect to the tangent to the isoline, namely by \(\langle \nabla g, J \nabla f \rangle\). If the curvature in the direction of the isoline \(k_{J \nabla f}\) is large compared to the stitch width \(w\), a uniform sampling rate is inadequate and does not result in a similar geometry. 

We therefore adjust the sampling rate in these regions by setting:
\[
\langle \nabla g, J \nabla f \rangle = h(k_{J \nabla f}),
\]
where \(h(x) = \frac{\tanh(-x / \alpha)}{2} + 1\), and \(\alpha = 10\). Figure 6 shows an example of such a model. In the central region, the model has negative mean and Gaussian curvature (a, b). Using a uniform sampling rate does not fully reconstruct the negative curvature (d, f), whereas a curvature-adapted sampling does (c, e).

---

### 7.2 Branching

If for a given seed \(s\), the geodesic function \(f(x) = d(s, x)\) has multiply-connected isolines, the graph \(G\) cannot cover the model. In these cases, the model is automatically decomposed into multiple segments, each of which can be crocheted using the approach described in Sections 4 and 5. The segments are attached by a method called **"join-as-you-go"** [Bennett 2020], meaning each segment is crocheted onto the last row of the previous segment, and therefore no additional sewing is required. Furthermore, the segment boundaries are not visible in the crocheted object. The more common method, in contrast, involves sewing together closed segments, which requires accuracy to achieve geometric details such as symmetry. 

#### 7.2.1 Mesh to Graph
Given a seed \(s\), let \(f(x) = d(s, x)\). Let \(\Pi = (\sigma_1, \ldots, \sigma_m)\) be the saddle points of \(f\), sorted by \(f_i = f(\sigma_i)\). Namely, \(f_1 \leq f_2 \leq \ldots \leq f_m\). For each \(\sigma_i\) in order, we compute the isoline of \(f_i\), denoted by \(\gamma_i\), and slice a new segment for each connected component of \(\gamma_i\). Meaning, the segments are obtained by slicing along the isolines of the saddle points, ordered by increasing geodesic distance to the seed. Figure 7(a) shows the isolines of \(f\) for the seed point \(s\) marked in red, as well as the saddles (cyan), and maxima (blue).

In addition to the segmentation, we generate a directed graph \(G_\sigma\), whose vertices are the segments, and where an edge \((s, t)\) exists if the segments \(M_s\) and \(M_t\) share a boundary, and \(f\) values on \(M_s\) are smaller than \(f\) values on \(M_t\). The crocheting order of the segments is determined by a topological sort of \(G_\sigma\). Figure 7(b) shows the resulting segments. Very thin segments might not be sampled (marked in black in Figure 7(b) and other figures) and are skipped and not crocheted.

**Geometry.** Any resulting segment \(M_l\) is either a half-sphere or a cylinder and thus can be covered by a crochet graph \(G_l\). While \(f\) is computed for the whole model before segmentation, \(g\) is computed for each segment separately. The cut for \(g\) is made from the maxima of \(f\) in the segment to the closest point on the segment’s boundary. If \(f\) attains its maximum on one of the boundaries of the segment (there are at most two boundaries), then the cut is computed to the closest point on the other boundary.

Figure 8 shows an example of the location of the cut, where we show the front and the back of the model. We show the location of the cut in brown on the segmented model (a), as well as the location of the cut in the crocheted model (b). To mark the location of the cut during crocheting, a piece of yarn was used to mark the beginning/end of the row.

**Connectivity.** For every two segments \(M_s, M_t\) which share an edge \((s, t)\) in \(G_\sigma\), we add an additional condition that the last row of \(M_s\) is coupled to the first row of \(M_t\). Figure 7(c) shows the crochet graph for all the segments of the Homer model. Finally, (e) shows the crocheted model, where each segment was crocheted with a different color for better visualization.

#### 7.2.2 Graph to Instructions
The same simulation of the crochet operations is applied to the first row of a new segment, but the sequence of stitches that is used as its previous row is no longer a full row. Instead, the last rows of all attached segments are arranged and filtered to include only vertices that have a connecting edge to the new segment’s row, constituting a “joint” previous row. The transducer then takes stock of when its consumption of the previous row skips stitches, splits segments, skips segments, or spans multiple parent segments, and includes this information in the row instructions.

---

**Figure 2 Description:**

This figure illustrates the **Crochet Graph** of a sphere, the corresponding crochet instructions (pattern), and the final crocheted object. It is divided into three parts:

1. **Left: The Crochet Graph**  
   The crochet graph of a sphere is shown. The graph is composed of:
   - **Row edges (\(R\))**: Highlighted in red, these edges connect consecutive stitches along each row of the sphere.
   - **Column edges (\(C\))**: Highlighted in blue, these edges connect stitches between adjacent rows.
   The sphere's graph is structured to ensure a grid-like pattern with approximately equal spacing between the stitches, following the curvature of the sphere.

2. **Middle: Crochet Instructions (Pattern)**  
   The crochet instructions to create the sphere are detailed row by row. Each line specifies the sequence of crochet operations, where:
   - `sc`: Single crochet.
   - `inc`: Increase stitch, which adds stitches to expand the row.
   - `dec`: Decrease stitch, which reduces stitches to taper the row.

   The instructions are as follows:
   - **Row 1**: 7sc in a ring [7 stitches total].
   - **Row 2**: (inc, sc)*2, 3inc [12 stitches total].
   - **Row 3**: inc, sc, (inc, 2sc)*2, (inc, sc)*2 [17 stitches total].
   - **Row 4**: 2sc, inc, 3sc, inc, 2sc, (2sc, inc, sc)*2 [21 stitches total].
   - **Row 5**: 4sc, inc, 10sc, inc, 5sc [23 stitches total].
   - **Row 6**: 12sc, inc, 10sc [24 stitches total].
   - **Row 7**: 24sc [24 stitches total].
   - **Row 8**: 3sc, dec, 6sc, dec, 7sc, dec, 2sc [21 stitches total].
   - **Row 9**: sc, dec, 3sc, dec, sc, (3sc, dec)*2, 2sc [17 stitches total].
   - **Row 10**: sc, dec, (sc, dec, sc)*2, (sc, dec)*2 [12 stitches total].
   - **Row 11**: dec, (dec, sc, dec)*2 [7 stitches total].
   - **Fasten off**: The process of finishing the crocheted sphere.

3. **Right: Final Crocheted Sphere**  
   A photograph of the crocheted sphere created using the instructions provided in the middle section. The sphere is composed of rows of single crochet stitches, with increases and decreases applied as specified to shape the sphere accurately. The fabric appears uniform and tightly woven, typical of amigurumi-style crochet.

This figure effectively demonstrates the transition from a mathematical representation (graph), to step-by-step instructions, to the tangible, crocheted outcome.

---

**Figure 3 Description:**

This figure demonstrates the structure and anatomy of crochet stitches using a visualized crocheted patch and a close-up analysis. The figure is divided into three panels:

1. **Left Panel: Crochet Stitches with Increase (\(inc\)) Marked**  
   A crocheted patch is shown with an **increase stitch (\(inc\))** illustrated. Increase stitches are marked with contrasting yarn (beige) woven into the patch. The increase stitch creates additional stitches in the row, allowing the fabric to expand in width.

2. **Middle Panel: Crochet Stitches with Decrease (\(dec\)) Marked**  
   Another crocheted patch is displayed, this time highlighting a **decrease stitch (\(dec\))**, again marked with contrasting beige yarn. Decrease stitches remove stitches from the row, causing the fabric to taper and narrow.

3. **Right Panel: Anatomy of a Crochet Stitch**  
   A close-up view of the crochet stitch structure is provided. The components of a single stitch are visually delineated with color-coded rectangles:
   - **Top/Bottom of the Stitch**: Highlighted with a red rectangle, this represents the points where stitches connect at their tops and bases. 
   - **Stem of the Stitch**: Highlighted with a green rectangle, this represents the vertical strand of yarn connecting the top of one stitch to the base of the next.

These visualizations and annotations clarify the functional differences between increase and decrease stitches and provide a detailed breakdown of the physical components of a crochet stitch.

---

**Figure 4 Description:**

This figure demonstrates the steps involved in generating a crochet graph \( G \) from the parameterization of a 3D shape. The figure is divided into four panels:

1. **Panel (a): The Isolines of \( f \) on the Shape \( M_C \)**  
   This panel shows the shape \( M_C \), a toroidal object, with the isolines of the function \( f \) visualized as evenly spaced, multicolored bands encircling the shape. Key annotations include:
   - **Seed Point \( s \)**: Marked in red, this point serves as the starting vertex for the parameterization process.
   - **Maxima of \( f \)**: Indicated in blue, representing the farthest point from the seed along the geodesic distance.
   - **Cut (black line)**: A geodesic cut along the shape, used to convert the toroidal topology into a topology suitable for grid sampling.
   The isolines represent rows of the crochet graph and reflect the progression of the parameter \( f \) from the seed point to the maxima.

2. **Panel (b): The \( f, g \) Parameterization and Sampled Grid \( S \)**  
   This panel shows the parameterization of the shape \( M_C \) into a 2D domain using the functions \( f \) and \( g \). The sampled grid \( S \) is overlaid on the domain, represented as an array of pink dots. The axes represent the values of \( f \) and \( g \), which correspond to the row and column indices in the crochet graph. The uniform distribution of points in the grid reflects the structured sampling of the parameterized surface.

3. **Panel (c): The Pushed Forward Points \( X_G \)**  
   The points from the sampled grid \( S \) are pushed forward onto the 3D shape \( M_C \) to generate the embedded vertices \( X_G \) of the crochet graph. The vertices are shown as a uniform distribution of pink dots covering the surface of the torus, maintaining the grid structure while conforming to the 3D geometry of \( M_C \).

4. **Panel (d): The Output Crochet Graph \( G \)**  
   The final crochet graph \( G \) is overlaid on the 3D shape \( M_C \). The graph consists of:
   - **Row Edges (\( R \))**: Highlighted in red, these edges connect consecutive stitches along each row.
   - **Column Edges (\( C \))**: Highlighted in blue, these edges connect stitches between adjacent rows.
   The grid-like structure of \( G \) approximates the geometry of \( M_C \), forming a regular crochet pattern adapted to the shape.

This figure illustrates the entire process of parameterization, sampling, and crochet graph generation, demonstrating how a complex 3D shape can be represented as a structured crochet graph for fabrication.

---

**Figure 5 Description:**

This figure illustrates the state of the transducer before and after producing a crochet stitch. It demonstrates different scenarios for the transitions between two rows of vertices, labeled \( S_i \) (current row) and \( S_{i+1} \) (next row). The figure is divided into four cases: `sc` (single crochet), `inc2` (increase by 2 stitches), `dec` (decrease stitch), and `illegal`. Each case shows the transducer state **before** (top row) and **after** (bottom row) processing. Key annotations include:

- **Green Circles**: Mark vertices at the heads of the rows being processed.
- **Purple Circles**: Represent unconsumed vertices in the rows.
- **Yellow Circles**: Indicate vertices that have been consumed during the stitching process.
- **Red Dotted Lines**: Denote the rows of the crochet graph.
- **Blue Lines**: Represent the column edges, connecting stitches between consecutive rows.

### Cases:

1. **`sc` (Single Crochet)**:
   - **Before**: One vertex from \( S_i \) is connected to one vertex in \( S_{i+1} \) with a single column edge (blue line). The vertices at the heads of the rows are highlighted in green.
   - **After**: Both vertices (one from \( S_i \) and one from \( S_{i+1} \)) are consumed and no longer at the head of the rows. The newly formed edge reflects a single crochet stitch.

2. **`inc2` (Increase by 2 Stitches)**:
   - **Before**: One vertex from \( S_i \) is connected to three vertices in \( S_{i+1} \), forming multiple column edges (blue lines).
   - **After**: The vertex in \( S_i \) is consumed, while two vertices in \( S_{i+1} \) remain unconsumed. This reflects an increase operation where two additional stitches are added to expand the row.

3. **`dec` (Decrease Stitch)**:
   - **Before**: Two vertices in \( S_i \) are connected to one vertex in \( S_{i+1} \), forming a single column edge (blue line).
   - **After**: The two vertices in \( S_i \) are consumed and merged into a single stitch in \( S_{i+1} \), representing a decrease operation that reduces the number of stitches in the row.

4. **`illegal` (Illegal Configuration)**:
   - **Before**: Two vertices in \( S_i \) are improperly connected to two vertices in \( S_{i+1} \), forming crossed column edges (blue lines) that violate the coupling constraints of the crochet graph.
   - **After**: This state is invalid as it results in an ambiguous or non-crochetable configuration. No legal transition exists from this state.

### Summary:
This figure illustrates the operation of the transducer, showing how vertices are processed to produce crochet stitches and highlighting legal and illegal configurations. Each transition ensures consistency with the coupling rules of the crochet graph to produce valid crochet instructions.

---

**Figure 6 Description:**

This figure demonstrates the effects of curvature-adapted and uniform sampling rates on a model with both negative mean curvature and negative Gaussian curvature in the same region. It consists of six panels:

1. **Panel (a): Negative Mean Curvature (\(H\)) Map**  
   This panel shows the mean curvature distribution over the 3D shape. The curvature values are visualized using a color gradient:
   - **Red**: Positive curvature.
   - **Blue**: Negative curvature.
   Areas with negative mean curvature are prominently visible, representing regions where the shape curves inward.

2. **Panel (b): Negative Gaussian Curvature (\(K\)) Map**  
   This panel illustrates the Gaussian curvature distribution over the same shape. The curvature values are visualized using a separate color gradient:
   - **Red**: Positive Gaussian curvature.
   - **Blue**: Negative Gaussian curvature.
   Regions with negative Gaussian curvature coincide with the areas where the surface exhibits a saddle-like shape.

3. **Panel (c): Curvature-Adapted Sampling (\(Y_G\))**  
   This panel shows the crochet graph generated using a curvature-adapted sampling rate. The crochet graph adheres more closely to the shape's geometry, with rows and columns of the graph following the curvature of the surface.

4. **Panel (d): Uniform Sampling (\(Y_G\))**  
   This panel displays the crochet graph generated using a uniform sampling rate. The uniform graph deviates from the input geometry in regions with high curvature, leading to inaccuracies in capturing the shape's features.

5. **Panel (e): Crocheted Object (Curvature-Adapted Sampling)**  
   This panel shows the crocheted object created using the curvature-adapted sampling rate. The shape closely resembles the input model, demonstrating the effectiveness of adapting the sampling rate based on curvature.

6. **Panel (f): Crocheted Object (Uniform Sampling)**  
   This panel illustrates the crocheted object created using a uniform sampling rate. The shape does not accurately replicate the input model, especially in regions with complex curvature.

### Key Observations:
- The curvature-adapted sampling method produces both a crochet graph and a crocheted object that closely resemble the input model's geometry.
- The uniform sampling method fails to accurately capture the shape's details in regions with high curvature, resulting in deviations from the intended geometry.

This figure emphasizes the importance of using curvature-adapted sampling for accurately representing and fabricating models with complex curvatures.

---

**Figure 7 Description:**

This figure illustrates the segmentation, crochet graph generation, and fabrication of a humanoid model. It consists of five panels, showing the progression from parameterization to the final crocheted object.

1. **Panel (a): Isolines of \( f \) on the Model**  
   The figure displays a humanoid model with colored isolines of the geodesic distance function \( f \), starting from the **seed point** marked in red. The isolines are evenly spaced, creating a striped appearance on the surface of the model. These isolines represent rows in the crochet graph.

2. **Panel (b): Segmented Model**  
   The humanoid model is divided into distinct **segments** based on the topology of \( f \). Saddle points, indicated in cyan, and maxima, indicated in blue, define the boundaries between segments. Each segment is shown in a different color to visually distinguish the components:
   - The torso, arms, and legs are separated into individual segments to facilitate crocheting.

3. **Panel (c): Crochet Graph \( G \)**  
   The crochet graph \( G \), consisting of vertices and edges, is overlaid on the segmented model. The graph includes:
   - **Row edges (\( R \))**: Highlighted in red, connecting consecutive stitches within the same row.
   - **Column edges (\( C \))**: Highlighted in blue, connecting stitches between adjacent rows.
   The graph captures the geometry and segmentation of the model, providing a structured plan for crocheting.

4. **Panel (d): Embedded Crochet Graph \( Y_G \)**  
   This panel shows the 3D embedding of the crochet graph \( Y_G \) on the surface of the model. The graph is evenly distributed over the shape, with edges closely aligned to the model's geometry. The embedded graph provides a direct representation of how the crochet pattern will align with the surface.

5. **Panel (e): Final Crocheted Model**  
   The fabricated humanoid model is displayed, crocheted using yarn in different colors for each segment. The segmentation allows for precise construction, with each part (e.g., head, arms, torso, legs) crocheted separately and seamlessly joined. The colors correspond to the segments shown in Panel (b), providing a clear visualization of the transition from graph to physical object.

### Key Observations:
- The isolines of \( f \) and segmentation ensure that the crochet graph adheres to the model's geometry and topology.
- The color-coded segments facilitate clarity in both the graph and the physical crocheted object.
- The final crocheted model accurately replicates the input geometry while showcasing the effectiveness of segmentation in crochet fabrication.

This figure highlights the complete pipeline, from parameterization and segmentation to graph generation and crochet fabrication.
---
**Figure 8 Description:**

This figure illustrates the placement of cuts on a segmented humanoid model for crocheting and the corresponding physical implementation. It is divided into two panels:

1. **Panel (a): Cut Location on the Segmented Model**  
   The humanoid model is displayed with distinct color-coded segments, each representing a different part of the body (e.g., head, arms, legs, torso). The location of the **cut** is highlighted in **brown** and runs along the back of the head and torso. Key details:
   - **Cut Line**: The cut separates the model to create an open topology for parameterization and crocheting. 
   - The cut starts at the top of the head near the **seed point** (marked in red) and continues downward along the center of the back, dividing the segments for crocheting.
   - Segments are color-coded to differentiate regions, such as the blue head, green torso, and yellow arms.

2. **Panel (b): Crocheted Model with Cut Marked**  
   This panel displays the physical crocheted humanoid model. The location of the cut is marked using a piece of contrasting yarn woven into the fabric. Key features include:
   - The **cut mark** (brown yarn) is clearly visible on the crocheted model, running vertically along the back to align with the cut location in Panel (a).
   - The crocheted object faithfully reproduces the shape and segmentation of the model, including smooth transitions between segments.
   - The left and right images show different views of the same crocheted object, illustrating the consistent placement of the cut mark.

### Key Observations:
- The placement of the cut facilitates the creation of a seamless and accurate crocheted object by dividing the model into manageable segments.
- The use of a physical marker (brown yarn) on the crocheted object ensures that the alignment of the cut location is clearly visible.
- The segmented approach results in a well-constructed crocheted model that preserves the geometry and proportions of the original design.

This figure demonstrates the correspondence between the digital segmentation and physical fabrication stages in the pipeline.

---
