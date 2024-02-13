#include <itensor/all_mps.h>
#include "CBE_TDVP.h"

int main() {
        int Ns = 10;

        int MaxDim = 2000;
        double Cutoff = 1e-10;
        double orth_cutoff = 1e-1;
        double expand_cutoff = 1e-3;
        double dt = 0.025;
        int Nt = 200;

        itensor::Args args("Cutoff", Cutoff, "MaxDim", MaxDim, "CutoffOrth", orth_cutoff, "CutoffExpand", expand_cutoff);

        auto sites = itensor::SpinHalf(Ns);
        auto ampo = itensor::AutoMPO(sites);
        for (int i=1; i < Ns; i++) {
                ampo += 0.5, "S+", i, "S-", i+1;
                ampo += 0.5, "S-", i, "S+", i+1;
                ampo += "Sz", i, "Sz", i+1;
        }
        double J2 = 0.5;
        for (int i=1; i < Ns-1; i++) {
                ampo += 0.5*J2, "S+", i, "S-", i+2;
                ampo += 0.5*J2, "S-", i, "S+", i+2;
                ampo += J2, "Sz", i, "Sz", i+2;
        }
        auto H = itensor::toMPO(ampo);
        auto state = itensor::InitState(sites);
        for (int i=1; i<= Ns; i++) {
                state.set(i, i%2 == 1 ? "Up" : "Dn");
        }
        auto psi = itensor::MPS(state);

        double nrm = 1.0;
        int center = Ns / 2;
        for (int i = 0; i <= Nt; i++) {
                std::vector<double> sz_tmp(Ns);
                psi.position(center);
                auto ket = psi.A(center);
                auto bra = itensor::dag(itensor::prime(ket, "Site"));
                double exp = itensor::eltC(ket * bra * itensor::op(sites, "Sz", center)).real();
                double ene = itensor::innerC(psi, H, psi).real();

                std::cout << "Time:" << i*dt << ", <Sz_center>:" << exp << ", Norm:" << nrm << ", Energy:" << ene << ", Dim:" << itensor::maxLinkDim(psi) << std::endl;
                nrm *= CBE_TDVP::CBE_TDVP(psi, H, dt, args);
        }

        return 0;
}
