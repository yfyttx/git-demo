// Exact GF(8) check of the generic G2 inner-factor/outer-check Turbo decoder.
#define main cascade_simulator_main
#include "simulate.cpp"
#undef main

int main() {
    GF2p gf(3,{1,0,1,1});
    FieldBasis basis(gf);
    NB_BP_Decoder decoder(gf,50);
    pauli_channel=PauliChannel::zBiased(0.12,1.0);

    // One generic inner generator row and one outer parity check.  The
    // resulting factor graph is a tree, so BP marginals must equal enumeration.
    constexpr int n1=2, n2=3, k2=1, q=8;
    SparseMatrix<int,RowMajor> outer(1,n1);
    outer.insert(0,0)=1;
    outer.insert(0,1)=1;
    const std::array<int,n2> g2={1,2,4};
    SparseMatrix<int,RowMajor> hx(k2,n1*n2);
    for(int block=0;block<n1;++block)
        for(int c=0;c<n2;++c)
            hx.insert(0,block*n2+c)=g2[c];
    hx.makeCompressed();

    const std::vector<std::vector<int>> x_profiles={
        {}, {0,1,2,4,3,5}, {7,6,5,4,3,2}
    };
    double max_error=0.0;
    double max_channel_mapping_error=0.0;
    int cases=0;
    int worst_profile=-1,worst_syndrome=-1,worst_variable=-1,worst_symbol=-1;
    double worst_bp=0.0,worst_exact=0.0;
    int profile_index=0;
    for(const auto& x_profile:x_profiles) {
    for(int syndrome_value=0;syndrome_value<q;++syndrome_value) {
        const std::vector<int>* x_estimate=x_profile.empty()?nullptr:&x_profile;
        std::vector<SoftVec> channel_priors(n1*n2,SoftVec(q,0.0));
        for(int v=0;v<n1*n2;++v) for(int z=0;z<q;++z) {
            double probability=1.0;
            const int xbits=x_estimate?basis.toDual((*x_estimate)[v]):0;
            for(int bit=0;bit<gf.getPower();++bit) {
                const double p=x_estimate
                    ? pauli_channel.zGivenX((xbits>>bit)&1)
                    : pauli_channel.marginalZ();
                probability*=((z>>bit)&1)?p:1.0-p;
            }
            channel_priors[v][z]=probability;
        }
        std::vector<int> syndrome={syndrome_value};
        std::vector<SoftVec> posterior;
        const auto estimate=decodeHxCascadeTurbo(
            decoder,outer,hx,syndrome,0.12,nullptr,basis,gf,
            12,20,0.8,false,&posterior,false,nullptr,false,nullptr,
            false,false,false,nullptr,false,&channel_priors);
        (void)estimate; // APP hard decisions need not satisfy the syndrome.

        if(x_estimate) {
            std::vector<SoftVec> mapped_posterior;
            (void)decodeHxCascadeTurbo(decoder,outer,hx,syndrome,0.12,
                x_estimate,basis,gf,12,20,0.8,false,&mapped_posterior);
            for(int v=0;v<n1*n2;++v) for(int z=0;z<q;++z)
                max_channel_mapping_error=std::max(max_channel_mapping_error,
                    std::abs(posterior[v][z]-mapped_posterior[v][z]));
        }

        std::vector<SoftVec> exact(n1*n2,SoftVec(q,0.0));
        double normalizer=0.0;
        for(int a=0;a<q;++a) for(int b=0;b<q;++b) for(int c=0;c<q;++c)
        for(int d=0;d<q;++d) for(int e=0;e<q;++e) for(int f=0;f<q;++f) {
            const std::array<int,n1*n2> word={a,b,c,d,e,f};
            int calculated=0;
            double probability=1.0;
            for(int block=0;block<n1;++block) {
                int y=0;
                for(int col=0;col<n2;++col) {
                    const int symbol=word[block*n2+col];
                    y=gf.add(y,gf.mul(g2[col],symbol));
                    probability*=channel_priors[block*n2+col][symbol];
                }
                calculated=gf.add(calculated,y);
            }
            if(calculated!=syndrome_value) continue;
            normalizer+=probability;
            for(int v=0;v<n1*n2;++v) exact[v][word[v]]+=probability;
        }
        for(int v=0;v<n1*n2;++v) for(int symbol=0;symbol<q;++symbol) {
            const double expected=exact[v][symbol]/normalizer;
            const double error=std::abs(posterior[v][symbol]-expected);
            if(error>max_error) {
                max_error=error;
                worst_profile=profile_index;
                worst_syndrome=syndrome_value;
                worst_variable=v;
                worst_symbol=symbol;
                worst_bp=posterior[v][symbol];
                worst_exact=expected;
            }
        }
        ++cases;
    }
    ++profile_index;
    }
    std::cout << cases << " exact GF(8) syndrome/conditional-prior cases; "
              << "max generic-G2 marginal error=" << max_error
              << " max x->z prior-map error=" << max_channel_mapping_error
              << " at profile=" << worst_profile << " syndrome="
              << worst_syndrome << " variable=" << worst_variable
              << " symbol=" << worst_symbol << " bp=" << worst_bp
              << " exact=" << worst_exact << '\n';
    return max_error<1e-8 && max_channel_mapping_error<1e-12 ? 0 : 1;
}
