#pragma once
#include <itensor/all_mps.h>

namespace ShrewdSelection {
        inline itensor::ITensor ShrewdSelectionRightToLeft(int index,
                                                    const itensor::ITensor& A, const itensor::ITensor& Lambda, const itensor::ITensor& B,
                                                    const itensor::MPO &H, const itensor::LocalMPO &PH, itensor::Args& args) {
                itensor::Index leftLink = index == 2 ? itensor::Index() : itensor::findIndex(A, tinyformat::format("Link,l=%d", index-2));
                itensor::Index siteInd = itensor::findIndex(A, "Site");
                itensor::Index lambdaLeft = itensor::commonIndex(A, Lambda);
                double orth_cutoff = args.getReal("CutoffOrth", 1e-4);
                double exp_cutoff = args.getReal("CutoffExpand", 1e-6);

                itensor::ITensor Rtmp = Lambda * B;
                if (PH.R()) {
                        Rtmp *= PH.R();
                }
                Rtmp *= H(index);
                itensor::ITensor Rorth = itensor::noPrime(Rtmp, "Site");
                Rtmp *= itensor::prime(itensor::dag(B));
                Rtmp *= itensor::prime(B, "Link");
                Rorth -= Rtmp;

                auto[RU, RS, RV] = itensor::svd(Rorth, {lambdaLeft}, args);
                auto Ltmp = A * RU * RS;
                if (PH.L()) {
                        Ltmp *= PH.L();
                }
                Ltmp *= H(index-1);
                itensor::ITensor Lorth = itensor::noPrime(Ltmp, "Site");
                Ltmp *= itensor::prime(itensor::dag(A));
                Ltmp *= itensor::prime(A, "Link");
                Lorth -= Ltmp;
                Lorth.noPrime();
                auto[LU, LS, LV] = itensor::svd(Lorth, {leftLink, siteInd, itensor::leftLinkIndex(H, index)}, {"Cutoff", orth_cutoff});
                auto[hatU, hatS, hatV] = itensor::svd(LU*LS, {leftLink, siteInd}, {"Cutoff", 1e-14});
                itensor::ITensor Apr = hatU;
                itensor::Index AprRight = itensor::commonIndex(hatU, hatS);

                itensor::ITensor Lpr = A;
                if (PH.L()) {
                        Lpr *= PH.L();
                }
                Lpr *= H(index-1);
                Lpr *= itensor::dag(itensor::prime(Apr));

                itensor::ITensor Ctmp = Lpr * Lambda * B;
                Ctmp *= H(index);
                if (PH.R()) {
                        Ctmp *= PH.R();
                }
                itensor::ITensor Corth = itensor::noPrime(Ctmp, "Site");
                Ctmp *= itensor::prime(itensor::dag(B));
                Ctmp *= itensor::prime(B, "Link");
                Corth -= Ctmp;

                auto [tildeU, tildeS, tildeV] = itensor::svd(Corth, {itensor::prime(AprRight)}, {"Cutoff", exp_cutoff, "LeftTags", tinyformat::format("Link,l=%d", index-1)});

                return Apr * itensor::noPrime(tildeU);
        }

        inline itensor::ITensor ShrewdSelectionLeftToRight(int index,
                                                    const itensor::ITensor& A, const itensor::ITensor& Lambda, const itensor::ITensor& B,
                                                    const itensor::MPO &H, const itensor::LocalMPO &PH, itensor::Args& args) {
                itensor::Index leftLink = index == 1 ? itensor::Index() : itensor::findIndex(A, tinyformat::format("Link,l=%d", index-1));
                itensor::Index wLink = itensor::rightLinkIndex(H, index);
                itensor::Index siteInd = itensor::findIndex(A, "Site");
                itensor::Index lambdaRight = itensor::commonIndex(B, Lambda);
                double orth_cutoff = args.getReal("CutoffOrth", 1e-4);
                double exp_cutoff = args.getReal("CutoffExpand", 1e-6);

                itensor::ITensor Ltmp = Lambda * A;
                if (PH.L()) {
                        Ltmp *= PH.L();
                }
                Ltmp *= H(index);
                itensor ::ITensor Lorth = itensor::noPrime(Ltmp, "Site");
                Ltmp *= itensor::prime(itensor::dag(A));
                Ltmp *= itensor::prime(A, "Link");
                Lorth -= Ltmp;
                Lorth.noPrime();

                auto[LU, LS, LV] = itensor::svd(Lorth, {leftLink, siteInd, wLink}, args);
                auto Rtmp = LS * LV * B;
                if (PH.R()) {
                        Rtmp *= PH.R();
                }
                Rtmp *= H(index+1);
                itensor::ITensor Rorth = itensor::noPrime(Rtmp, "Site");
                Rtmp *= itensor::prime(itensor::dag(B));
                Rtmp *= itensor::prime(B, "Link");
                Rorth -= Rtmp;
                Rorth.noPrime();
                auto[RU, RS, RV] = itensor::svd(Rorth, {itensor::commonIndex(LU, LS)}, {"Cutoff", orth_cutoff});
                auto[hatU, hatS, hatV] = itensor::svd(RS*RV, {itensor::commonIndex(RU, RS), wLink}, {"Cutoff", 1e-14});
                itensor::ITensor Bpr = hatV;

                itensor::ITensor Rpr = B;
                if (PH.R()) {
                        Rpr *= PH.R();
                }
                Rpr *= H(index+1);
                Rpr *= itensor::dag(itensor::prime(Bpr));

                itensor::ITensor Ctmp = Rpr * Lambda * A;
                Ctmp *= H(index);
                if (PH.L()) {
                        Ctmp *= PH.L();
                }
                itensor::ITensor Corth = itensor::noPrime(Ctmp, "Site");
                Ctmp *= itensor::prime(itensor::dag(A));
                Ctmp *= itensor::prime(A, "Link");
                Corth -= Ctmp;
                Corth.noPrime();

                auto [tildeU, tildeS, tildeV] = itensor::svd(Corth, {leftLink, siteInd}, {"Cutoff", exp_cutoff, "RightTags", tinyformat::format("Link,l=%d", index)});

                return Bpr * itensor::noPrime(tildeV);
        }
}
