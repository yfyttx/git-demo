// Exact GF(8) enumeration for the constrained physical-symbol marginals.
#define main cascade_simulator_main
#include "simulate.cpp"
#undef main

int main() {
    GF2p gf(3,{1,0,1,1});
    FieldBasis basis(gf);
    NB_BP_Decoder decoder(gf,50);
    pauli_channel=PauliChannel::zBiased(0.12,1.0);
    SparseMatrix<int,RowMajor> outer(1,1); outer.insert(0,0)=1;
    double maximum_error=0;
    for(bool weighted:{false,true}) for(int forced:{-1,0,1})
        for(bool conditional:{false,true}) {
        std::vector<int> x={1,2,4};
        std::vector<int> scale={1,weighted?gf.exp(1):1,weighted?gf.exp(2):1};
        SparseMatrix<int,RowMajor> hx(2,3);
        hx.insert(0,0)=scale[0]; hx.insert(0,1)=scale[1];
        hx.insert(1,1)=scale[1]; hx.insert(1,2)=scale[2];
        std::vector<int> syn={3,1};
        std::vector<std::pair<int,int>> constraints(3,{0,0});
        if(forced>=0) constraints[0]={1,forced};
        std::vector<SoftVec> post;
        auto decoded=decodeHxCascadeTurbo(decoder,outer,hx,syn,0.12,conditional ? &x : nullptr,
            basis,gf,1,1,0.2,false,&post,false,nullptr,false,nullptr,false,false,
            false,&constraints,true);
        if(!decoder.check_z_syndrome(hx,decoded,syn))
            throw std::runtime_error("Invalid reconstructed syndrome");
        std::vector<SoftVec> exact(3,SoftVec(8,0));
        double total=0;
        for(int a=0;a<8;++a) for(int b=0;b<8;++b) for(int c=0;c<8;++c) {
            std::vector<int> z={a,b,c};
            if(!decoder.check_z_syndrome(hx,z,syn)) continue;
            double probability=1;
            for(int v=0;v<3;++v) for(int bit=0;bit<3;++bit) {
                const double pz=conditional
                    ? pauli_channel.zGivenX((basis.toDual(x[v])>>bit)&1) : 0.08;
                probability*=((z[v]>>bit)&1)?pz:1-pz;
            }
            if(forced>=0 && (a&1)!=forced) probability*=std::exp(-25.0);
            total+=probability;
            for(int v=0;v<3;++v) exact[v][z[v]]+=probability;
        }
        for(int v=0;v<3;++v) for(int z=0;z<8;++z)
            maximum_error=std::max(maximum_error,std::abs(post[v][z]-exact[v][z]/total));
    }
    std::cout << "12 exact GF(8) cases; maximum marginal error=" << maximum_error << '\n';
    if(maximum_error>1e-8) return 1;
}
