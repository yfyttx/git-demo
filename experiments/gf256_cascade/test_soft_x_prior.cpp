#define main gf256_cascade_simulator_main
#include "simulate.cpp"
#undef main

int main() {
    GF2p gf(3,{1,0,1,1});
    FieldBasis basis(gf);
    NB_BP_Decoder decoder(gf,50);
    pauli_channel=PauliChannel::zBiased(0.21,2.5);

    // Exact GF(8) check for the SI local solve: the dual row [1,1,0]
    // lies in ker(H), so its first two variables can be inactivated. The
    // remaining degree-one check fixes x2; enumerate all local solutions
    // and verify that SI selects the posterior-maximizing one.
    {
        SparseMatrix<int,RowMajor> si_h(2,3),si_dual(1,3);
        si_h.insert(0,0)=1; si_h.insert(0,1)=1;
        si_h.insert(1,2)=1;
        si_dual.insert(0,0)=1; si_dual.insert(0,1)=1;
        si_h.makeCompressed(); si_dual.makeCompressed();
        const std::vector<int> si_truth={2,3,4};
        std::vector<int> si_syndrome(2,0);
        for(int r=0;r<si_h.rows();++r)
            for(SparseMatrix<int,RowMajor>::InnerIterator it(si_h,r);it;++it)
                si_syndrome[r]=gf.add(si_syndrome[r],
                                      gf.mul(it.value(),si_truth[it.col()]));
        std::vector<SoftVec> si_priors(3,SoftVec(gf.getSize(),
                                                   1.0/gf.getSize()));
        std::vector<SoftVec> si_posterior(3,SoftVec(gf.getSize(),0.0));
        for(int x=0;x<gf.getSize();++x) {
            si_posterior[0][x]=x==2?0.90:0.10/(gf.getSize()-1);
            si_posterior[1][x]=x==3?0.95:0.05/(gf.getSize()-1);
            si_posterior[2][x]=1.0/gf.getSize();
        }
        const auto si=decodeOuterWithStabilizerInactivation(
            decoder,si_h,si_dual,si_syndrome,si_priors,si_posterior,8,1,gf);
        if(!si.syndrome_ok || !decoder.check_z_syndrome(
                si_h,si.estimate,si_syndrome))
            throw std::runtime_error("GF(8) SI failed to recover the syndrome");
        double best_score=-std::numeric_limits<double>::infinity();
        std::vector<int> best_pair;
        for(int x0=0;x0<gf.getSize();++x0)
            for(int x1=0;x1<gf.getSize();++x1)
                if(gf.add(x0,x1)==si_syndrome[0]) {
                    const double score=std::log(si_posterior[0][x0])+
                                       std::log(si_posterior[1][x1]);
                    if(score>best_score) {
                        best_score=score;
                        best_pair={x0,x1};
                    }
                }
        if(si.estimate[0]!=best_pair[0] || si.estimate[1]!=best_pair[1] ||
           si.estimate[2]!=si_truth[2])
            throw std::runtime_error("GF(8) SI disagrees with exhaustive local MAP");
        std::cout << "GF(8) SI exact local solve passed; attempts="
                  << si.attempts << '\n';
    }

    SparseMatrix<int,RowMajor> h(2,3);
    h.insert(0,0)=1;
    h.insert(0,1)=gf.exp(1);
    h.insert(1,1)=gf.exp(2);
    h.insert(1,2)=1;
    h.makeCompressed();
    const std::vector<int> truth={1,4,6};
    std::vector<int> syndrome(2,0);
    for(int r=0;r<2;++r)
        for(SparseMatrix<int,RowMajor>::InnerIterator it(h,r);it;++it)
            syndrome[r]=gf.add(syndrome[r],gf.mul(it.value(),truth[it.col()]));

    const int q=gf.getSize();
    const double px=pauli_channel.marginalX();
    SoftVec channel_prior(q,0.0);
    for(int x=0;x<q;++x) {
        const int bits=basis.toDual(x);
        double probability=1.0;
        for(int bit=0;bit<gf.getPower();++bit)
            probability*=((bits>>bit)&1)?px:1.0-px;
        channel_prior[x]=probability;
    }
    std::vector<SoftVec> priors(3,channel_prior);
    const auto bp=runOuterSoftBP(decoder,h,syndrome,priors,10,gf);
    if(!bp.valid) throw std::runtime_error("Tree GF(8) BP did not satisfy syndrome");

    std::vector<SoftVec> exact_x(3,SoftVec(q,0.0));
    double partition=0.0;
    for(int x0=0;x0<q;++x0) for(int x1=0;x1<q;++x1) for(int x2=0;x2<q;++x2) {
        const std::vector<int> x={x0,x1,x2};
        bool valid=true;
        for(int r=0;r<2;++r) {
            int value=0;
            for(SparseMatrix<int,RowMajor>::InnerIterator it(h,r);it;++it)
                value=gf.add(value,gf.mul(it.value(),x[it.col()]));
            valid &= value==syndrome[r];
        }
        if(!valid) continue;
        const double weight=channel_prior[x0]*channel_prior[x1]*channel_prior[x2];
        partition+=weight;
        for(int v=0;v<3;++v) exact_x[v][x[v]]+=weight;
    }
    double max_x_error=0.0;
    for(int v=0;v<3;++v) for(int x=0;x<q;++x)
        max_x_error=std::max(max_x_error,
            std::abs(bp.posterior[v][x]-exact_x[v][x]/partition));

    const auto transformed=softXPosteriorToZPrior(
        bp.posterior,basis,gf,pauli_channel);
    double max_z_error=0.0;
    for(int v=0;v<3;++v) for(int z=0;z<q;++z) {
        double expected=0.0;
        for(int x=0;x<q;++x) {
            const int xbits=basis.toDual(x);
            double conditional=1.0;
            for(int bit=0;bit<gf.getPower();++bit) {
                const double pz=pauli_channel.zGivenX((xbits>>bit)&1);
                conditional*=((z>>bit)&1)?pz:1.0-pz;
            }
            expected+=bp.posterior[v][x]*conditional;
        }
        max_z_error=std::max(max_z_error,
            std::abs(transformed[v][z]-expected));
    }

    // Verify the Turbo decoder's scaled internal Z-symbol prior against the
    // physical polynomial-basis Z bits, conditioned on the trace-dual X bits.
    // The Turbo variable is scale * physical_z, so its likelihood must be
    // indexed by physical_z = internal_z / scale.
    {
        const std::vector<int> x_symbols={1,5,6};
        const std::vector<int> scales={1,gf.exp(1),gf.exp(2)};
        double max_scaled_prior_error=0.0;
        for(int v=0;v<3;++v) {
            SoftVec turbo_prior(q,0.0), physical_prior_reindexed(q,0.0);
            const int xbits=basis.toDual(x_symbols[v]);
            for(int internal_z=0;internal_z<q;++internal_z) {
                const int physical_z=gf.mul(internal_z,gf.inv(scales[v]));
                double probability=1.0;
                for(int bit=0;bit<gf.getPower();++bit) {
                    const double pz=pauli_channel.zGivenX((xbits>>bit)&1);
                    probability*=((physical_z>>bit)&1)?pz:1.0-pz;
                }
                turbo_prior[internal_z]=probability;
            }
            for(int physical_z=0;physical_z<q;++physical_z) {
                const int internal_z=gf.mul(scales[v],physical_z);
                double probability=1.0;
                for(int bit=0;bit<gf.getPower();++bit) {
                    const double pz=pauli_channel.zGivenX((xbits>>bit)&1);
                    probability*=((physical_z>>bit)&1)?pz:1.0-pz;
                }
                physical_prior_reindexed[internal_z]=probability;
            }
            normalizeLocal(turbo_prior);
            normalizeLocal(physical_prior_reindexed);
            for(int z=0;z<q;++z)
                max_scaled_prior_error=std::max(max_scaled_prior_error,
                    std::abs(turbo_prior[z]-physical_prior_reindexed[z]));
        }
        std::cout << "GF(8) scaled conditional-Z prior max error="
                  << max_scaled_prior_error << '\n';
        if(max_scaled_prior_error>1e-14)
            throw std::runtime_error("Scaled Z prior disagrees with physical-bit prior");
    }

    // Exact tree check for the structured selector + per-block H2 posterior.
    SparseMatrix<int,RowMajor> outer(1,2);
    outer.insert(0,0)=1;
    outer.makeCompressed();
    SparseMatrix<int,RowMajor> hz_cascade(4,6);
    hz_cascade.insert(0,0)=1;
    hz_cascade.insert(1,1)=1;
    hz_cascade.insert(2,0)=1; hz_cascade.insert(2,1)=1; hz_cascade.insert(2,2)=1;
    hz_cascade.insert(3,3)=1; hz_cascade.insert(3,4)=1; hz_cascade.insert(3,5)=1;
    hz_cascade.makeCompressed();
    const std::vector<int> x_truth={1,4,6,2,3,5};
    std::vector<int> x_syn(hz_cascade.rows(),0);
    for(int r=0;r<hz_cascade.rows();++r)
        for(SparseMatrix<int,RowMajor>::InnerIterator it(hz_cascade,r);it;++it)
            x_syn[r]=gf.add(x_syn[r],gf.mul(it.value(),x_truth[it.col()]));
    SoftVec bit_channel(q,0.0);
    for(int x=0;x<q;++x) {
        const int bits=basis.toDual(x);
        double prob=1.0;
        for(int b=0;b<gf.getPower();++b)
            prob*=((bits>>b)&1)?px:1.0-px;
        bit_channel[x]=prob;
    }
    std::vector<SoftVec> cascade_priors(6,bit_channel);
    const auto structured=runCascadeXSoftBP(
        decoder,outer,hz_cascade,x_syn,cascade_priors,1,2,2,gf);
    double max_structured_error=0.0;
    for(int block=0;block<2;++block) {
        std::vector<SoftVec> exact(3,SoftVec(q,0.0));

        for(int a=0;a<q;++a) for(int b=0;b<q;++b) for(int c=0;c<q;++c) {
            const std::vector<int> word={a,b,c};
            if(block==0 && (a!=x_syn[0] || b!=x_syn[1])) continue;
            if(gf.add(gf.add(a,b),c)!=x_syn[2+block]) continue;
            const double weight=bit_channel[a]*bit_channel[b]*bit_channel[c];

            exact[0][a]+=weight; exact[1][b]+=weight; exact[2][c]+=weight;
        }
        for(auto& marginal:exact) normalizeLocal(marginal);
        for(int local=0;local<3;++local)
            for(int value=0;value<q;++value)
                max_structured_error=std::max(max_structured_error,
                    std::abs(structured.posterior[block*3+local][value]-exact[local][value]));
    }

    std::cout << "GF(8) tree posterior max error=" << max_x_error
              << "; Pauli X-to-Z transform max error=" << max_z_error
              << "; structured cascade posterior max error=" << max_structured_error << '\n';
    if(max_x_error>1e-9 || max_z_error>1e-12 || max_structured_error>1e-9) return 1;
    return 0;
}
