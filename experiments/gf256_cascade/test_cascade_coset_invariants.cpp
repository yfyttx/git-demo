#define main gf256_simulator_main
#include "simulate.cpp"
#undef main

static std::vector<std::vector<int>> rightNullspaceBasis(
    const SparseMatrix<int,RowMajor>& matrix, GF2p& gf) {
    std::vector<std::vector<int>> a(matrix.rows(),
                                    std::vector<int>(matrix.cols(),0));
    for (int r=0;r<matrix.rows();++r)
        for (SparseMatrix<int,RowMajor>::InnerIterator it(matrix,r);it;++it)
            a[r][it.col()]=it.value();
    std::vector<int> pivots;
    int rank=0;
    for (int col=0;col<matrix.cols() && rank<matrix.rows();++col) {
        int pivot=rank;
        while (pivot<matrix.rows() && a[pivot][col]==0) ++pivot;
        if (pivot==matrix.rows()) continue;
        std::swap(a[rank],a[pivot]);
        const int inverse=gf.inv(a[rank][col]);
        for (int j=col;j<matrix.cols();++j)
            a[rank][j]=gf.mul(a[rank][j],inverse);
        for (int r=0;r<matrix.rows();++r) {
            if (r==rank || a[r][col]==0) continue;
            const int factor=a[r][col];
            for (int j=col;j<matrix.cols();++j)
                a[r][j]=gf.add(a[r][j],gf.mul(factor,a[rank][j]));
        }
        pivots.push_back(col);
        ++rank;
    }
    std::vector<bool> is_pivot(matrix.cols(),false);
    for (int pivot:pivots) is_pivot[pivot]=true;
    std::vector<std::vector<int>> basis;
    for (int free_col=0;free_col<matrix.cols();++free_col) {
        if (is_pivot[free_col]) continue;
        std::vector<int> v(matrix.cols(),0);
        v[free_col]=1;
        for (int r=0;r<rank;++r) v[pivots[r]]=a[r][free_col];
        basis.push_back(std::move(v));
    }
    return basis;
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const std::string dir = argv[1];
    const std::vector<int> poly{1,0,0,0,1,1,1,0,1};
    GF2p gf(8, poly);
    auto hx = loadMatrix(dir + "/HX_Cascade_GF256.txt", gf);
    auto hz = loadMatrix(dir + "/HZ_Cascade_GF256.txt", gf);
    auto outer = loadMatrix(dir + "/HX_TIT_Outer_For_HZ_GF256.txt", gf);
    auto hx_outer = loadMatrix(dir + "/HZ_TIT_Outer_For_HX_GF256.txt", gf);
    SparseMatrix<int,RowMajor> h2;
    SparseMatrix<int,RowMajor> g2;
    if (std::ifstream(dir + "/H2_Inner_GF256.txt"))
        h2 = loadMatrix(dir + "/H2_Inner_GF256.txt", gf);
    if (std::ifstream(dir + "/G2_Inner_GF256.txt"))
        g2 = loadMatrix(dir + "/G2_Inner_GF256.txt", gf);
    if (hx.cols() != hz.cols() || hx.cols() % 12 != 0 || outer.rows() == 0 ||
        hx.rows() % outer.rows() != 0) return 3;
    const int n2=12, k2=hx.rows()/outer.rows();
    std::vector<int> scales(n2,1);
    if (h2.rows()==1) {
        for (int j=0;j<n2;++j) {
            const int v=h2.coeff(0,j);
            if (!v) return 4;
            scales[j]=gf.inv(v);
        }
    }
    FieldRowSpace outer_space(outer,gf,h2.rows() ? &h2 : nullptr,
                              g2.rows() ? &g2 : nullptr);
    FieldRowSpace x_outer_space(hx_outer,gf,nullptr,g2.rows() ? &g2 : nullptr);
    int rejected_x_rows=0;
    if (!isCascadeXStabilizer(std::vector<int>(hx.cols(),0),
                             x_outer_space,gf,n2,k2,scales)) ++rejected_x_rows;
    for (int r=0;r<hx.rows();++r) {
        std::vector<int> row(hx.cols(),0);
        for (SparseMatrix<int,RowMajor>::InnerIterator it(hx,r);it;++it)
            row[it.col()]=it.value();
        if (!isCascadeXStabilizer(row,x_outer_space,gf,n2,k2,scales))
            ++rejected_x_rows;
    }
    FieldRowSpace hz_row_space(hz,gf);
    std::vector<std::vector<std::pair<int,int>>> hx_by_col(hx.cols());
    for (int r=0;r<hx.rows();++r)
        for (SparseMatrix<int,RowMajor>::InnerIterator it(hx,r);it;++it)
            hx_by_col[it.col()].push_back({r,it.value()});
    std::vector<int> deterministic(hx.cols());
    for (int i=0;i<hx.cols();++i) deterministic[i]=(i*73+19)%256;
    const auto ref_key=cascadeZCosetKey(deterministic,outer_space,gf,n2,k2,scales);
    int nonorth_rows=0, nonstabilizer_rows=0, key_changes=0;
    for (int r=0;r<hz.rows();++r) {
        std::vector<int> row(hz.cols(),0), syndrome(hx.rows(),0);
        for (SparseMatrix<int,RowMajor>::InnerIterator it(hz,r);it;++it) {
            row[it.col()]=it.value();
            for (const auto& [hr,v] : hx_by_col[it.col()])
                syndrome[hr]=gf.add(syndrome[hr],gf.mul(v,it.value()));
        }
        if (std::any_of(syndrome.begin(),syndrome.end(),[](int x){return x!=0;})) ++nonorth_rows;
        if (!isCascadeZStabilizer(row,outer_space,gf,n2,k2,scales)) ++nonstabilizer_rows;
        auto shifted=deterministic;
        for (int j=0;j<hz.cols();++j) shifted[j]=gf.add(shifted[j],row[j]);
        if (cascadeZCosetKey(shifted,outer_space,gf,n2,k2,scales)!=ref_key) ++key_changes;
    }
    const auto kernel_basis=rightNullspaceBasis(hx_outer,gf);
    int found_consistent_logical_witness=0;
    std::vector<int> logical_witness(hx.cols(),0);
    for (const auto& coordinate:kernel_basis) {
        if (outer_space.contains(coordinate)) continue;
        std::fill(logical_witness.begin(),logical_witness.end(),0);
        for (int j=0;j<static_cast<int>(coordinate.size());++j)
            logical_witness[j*n2]=coordinate[j];
        std::vector<int> witness_syndrome(hx.rows(),0);
        for (int j=0;j<hx.cols();++j) if (logical_witness[j])
            for (const auto& [hr,v]:hx_by_col[j])
                witness_syndrome[hr]=gf.add(witness_syndrome[hr],
                    gf.mul(v,logical_witness[j]));
        const bool syndrome_zero=std::all_of(witness_syndrome.begin(),
            witness_syndrome.end(),[](int x){return x==0;});
        const bool helper_stabilizer=isCascadeZStabilizer(
            logical_witness,outer_space,gf,n2,k2,scales);
        const bool true_stabilizer=hz_row_space.contains(logical_witness);
        const auto zero_key=cascadeZCosetKey(
            std::vector<int>(hx.cols(),0),outer_space,gf,n2,k2,scales);
        const bool key_is_trivial=cascadeZCosetKey(
            logical_witness,outer_space,gf,n2,k2,scales)==zero_key;
        if (syndrome_zero && !helper_stabilizer && !true_stabilizer &&
            !key_is_trivial) {
            found_consistent_logical_witness=1;
            break;
        }
    }
    std::cout << "matrix_dir=" << dir << "\n"
              << "dimensions HX=" << hx.rows() << "x" << hx.cols()
              << " HZ=" << hz.rows() << "x" << hz.cols()
              << " outer=" << outer.rows() << "x" << outer.cols()
              << " H2=" << h2.rows() << "x" << h2.cols()
              << " inner_dimension=" << k2 << "\n"
              << "HZ_rows_nonorthogonal_to_HX=" << nonorth_rows << "/" << hz.rows() << "\n"
              << "HX_rows_and_zero_rejected_by_isCascadeXStabilizer="
              << rejected_x_rows << "/" << hx.rows()+1 << "\n"
              << "HZ_rows_rejected_by_isCascadeZStabilizer=" << nonstabilizer_rows << "/" << hz.rows() << "\n"
              << "HZ_rows_changing_coset_key=" << key_changes << "/" << hz.rows() << "\n"
              << "nontrivial_logical_witness_classified_consistently="
              << found_consistent_logical_witness << "\n";
    return (rejected_x_rows || nonorth_rows || nonstabilizer_rows || key_changes ||
            !found_consistent_logical_witness) ? 1 : 0;
}
