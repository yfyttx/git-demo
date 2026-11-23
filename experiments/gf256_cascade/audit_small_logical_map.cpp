// Exhaustive logical-MAP benchmark for a six-symbol GF(8) CSS cascade.
// H_X=[1,alpha] kron G2; H_Z=I_2 kron [1,1,1]. No X decoding:
// conditional profiles below are fixed side-information models, not FER trials.
#define main cascade_simulator_main
#include "simulate.cpp"
#undef main

static int class_key(const std::vector<int>& z) {
    int key=0;
    for(int b=0;b<2;++b) {
        key=8*key+(z[3*b]^z[3*b+2]);
        key=8*key+(z[3*b+1]^z[3*b+2]);
    }
    return key;
}
int main() {
    GF2p gf(3,{1,0,1,1}); FieldBasis basis(gf); NB_BP_Decoder decoder(gf,50);
    std::cout << "basis,p_D,profile,damping,rounds,logical_map_risk,error_map_logical_risk,exact_marginal_hard_risk,turbo_logical_risk,turbo_syndrome_failure_mass,suboptimal_syndromes,max_marginal_error,mean_symbol_total_variation\n";
    for(bool path : {false,true}) {
        SparseMatrix<int,RowMajor> outer(1,2),hx(2,6);
        outer.insert(0,0)=1; outer.insert(0,1)=2; outer.makeCompressed();
        for(int b=0;b<2;++b) for(int a=0;a<2;++a) {
            hx.insert(a,3*b+(path?a:0))=b?2:1;
            hx.insert(a,3*b+a+1)=b?2:1;
        }
        hx.makeCompressed();
        auto syn=[&](const std::vector<int>& z) {
            int s0=0,s1=0;
            for(int b=0;b<2;++b) {
                s0 ^= gf.mul(b?2:1,z[3*b]^z[3*b+1]);
                s1 ^= gf.mul(b?2:1,z[3*b+(path?1:0)]^z[3*b+2]);
            }
            return 8*s0+s1;
        };
        const std::vector<int> zero(6,0);
        for(int block=0;block<2;++block) {
            auto stabilizer=zero;
            for(int c=0;c<3;++c) stabilizer[3*block+c]=1;
            if(syn(stabilizer)!=syn(zero) ||
               class_key(stabilizer)!=class_key(zero))
                throw std::runtime_error("Toy CSS orthogonality/coset check failed");
        }
        for(double p : {0.04,0.07125,0.12}) for(int profile=0;profile<3;++profile) {
            pauli_channel=PauliChannel::zBiased(p,1.0); // pX=pY=pZ=p/3
            std::vector<SoftVec> prior(6,SoftVec(8));
            const int xbits[6]={1,2,0,4,0,3};
            for(int v=0;v<6;++v) for(int z=0;z<8;++z) {
                double mass=1;
                for(int bit=0;bit<3;++bit) {
                    double flip=profile==0 ? 2*p/3 :
                        pauli_channel.zGivenX(profile==1?0:((xbits[v]>>bit)&1));
                    mass*=((z>>bit)&1)?flip:1-flip;
                }
                prior[v][z]=mass;
            }
            std::vector<double> class_mass(4096,0),syndrome_mass(64,0),best_word_mass(64,-1);
            std::vector<int> class_syn(4096,-1),best_word_key(64,0);
            std::vector<std::vector<SoftVec>> marginal(64,std::vector<SoftVec>(6,SoftVec(8,0)));
            for(int packed=0;packed<(1<<18);++packed) {
                std::vector<int> z(6);int code=packed;double mass=1;
                for(int v=0;v<6;++v){z[v]=code&7;code>>=3;mass*=prior[v][z[v]];}
                int s=syn(z),k=class_key(z);
                if(class_syn[k]!=-1 && class_syn[k]!=s) throw std::runtime_error("Invalid stabilizer coset");
                class_syn[k]=s;class_mass[k]+=mass;syndrome_mass[s]+=mass;
                for(int v=0;v<6;++v)marginal[s][v][z[v]]+=mass;
                if(mass>best_word_mass[s]){best_word_mass[s]=mass;best_word_key[s]=k;}
            }
            std::vector<double> best_class_mass(64,0);
            for(int k=0;k<4096;++k)best_class_mass[class_syn[k]]=std::max(best_class_mass[class_syn[k]],class_mass[k]);
            double total=std::accumulate(syndrome_mass.begin(),syndrome_mass.end(),0.0);
            if(std::abs(total-1)>1e-10)throw std::runtime_error("Probability not normalized");
            double optimal=0,error_map=0,exact_hard=0;
            for(int s=0;s<64;++s){
                optimal+=best_class_mass[s];error_map+=class_mass[best_word_key[s]];
                std::vector<int> hard(6);
                for(int v=0;v<6;++v){
                    for(double& value:marginal[s][v])value/=syndrome_mass[s];
                    hard[v]=std::max_element(marginal[s][v].begin(),marginal[s][v].end())-marginal[s][v].begin();
                }
                if(syn(hard)==s)exact_hard+=class_mass[class_key(hard)];
            }
            for(auto settings : std::vector<std::pair<double,int>>{{0.2,20},{0.8,20},{0.2,60}}){
                double success=0,invalid=0,max_error=0,weighted_tv=0;int worse=0;
                for(int s=0;s<64;++s){
                    std::vector<SoftVec> post;
                    auto est=decodeHxCascadeTurbo(decoder,outer,hx,{s/8,s%8},p,nullptr,basis,gf,
                        settings.second,20,settings.first,false,&post,false,nullptr,false,nullptr,
                        false,false,false,nullptr,false,&prior);
                    double mass=0;
                    if(syn(est)==s)mass=class_mass[class_key(est)];else invalid+=syndrome_mass[s];
                    success+=mass;
                    if(best_class_mass[s]-mass>1e-12)++worse;
                    for(int v=0;v<6;++v)for(int z=0;z<8;++z)
                    {
                        const double difference=std::abs(post[v][z]-marginal[s][v][z]);
                        max_error=std::max(max_error,difference);
                        weighted_tv+=syndrome_mass[s]*difference/12.0;
                    }
                }
                if(success>optimal+1e-10 || error_map>optimal+1e-10)throw std::runtime_error("MAP bound violated");
                std::cout << (path?"path":"star") << ',' << p << ',' << profile << ','
                    << settings.first << ',' << settings.second << ',' << 1-optimal << ','
                    << 1-error_map << ',' << 1-exact_hard << ',' << 1-success << ','
                    << invalid << ',' << worse << ',' << max_error << ',' << weighted_tv << '\n';
            }
        }
    }
}
