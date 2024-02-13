#include <itensor/all_mps.h>
#include "CBE_DMRG.h"

int main() {
        int Ns = 100;
        auto sites = itensor::Boson(Ns, {"MaxOcc", 10});
        auto ini_state = itensor::InitState(sites, "1");
        auto psi = itensor::MPS(ini_state);

        double U = 2.0;

        auto ampo = itensor::AutoMPO(sites);
        for (int i = 1; i < Ns; i++) {
                ampo += "A", i , "Adag", i+1;
                ampo += "Adag", i , "A", i+1;
        }

        for (int i = 1; i <= Ns; i++) {
                ampo += U, "N*N", i;
                ampo += -U, "N", i;
        }
        double cutoff = 1e-10;
        int maxdim = 200;
        double cutoff_orth = 1e-1;
        double cutoff_expand = 1e-3;
        itensor::Args args(
                        "Cutoff", cutoff,
                        "MaxDim", maxdim,
                        "CutoffOrth", cutoff_orth,
                        "CutoffExpand", cutoff_expand
                        );
        auto H = itensor::toMPO(ampo);
        double ene = CBE_DMRG::CBE_DMRG(psi, H, args);
        double ene_prev = ene;
        while (true) {
                ene = CBE_DMRG::CBE_DMRG(psi, H, args);
                double rel_err = std::abs((ene_prev - ene) / ene);
                std::cout << "Max bond dim:" << itensor::maxLinkDim(psi) << std::endl;
                std::cout << "relative error:" << rel_err << std::endl;
                if (rel_err < cutoff) {
                        break;
                }
                ene_prev = ene;
        }
        std::cout << std::setprecision(16) << "obtained energy:" << ene << std::endl;
        return 0;
}
