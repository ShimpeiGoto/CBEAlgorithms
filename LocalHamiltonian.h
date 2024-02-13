#pragma once
#include <itensor/all_mps.h>

namespace LocalHamiltonian {
        class SingleSiteH {
                private:
                        const itensor::ITensor& LTensor_, RTensor_, W_;
                        size_t size_;
                public:
                        SingleSiteH(const itensor::ITensor& LTensor, const itensor::ITensor& RTensor, const itensor::ITensor& W) : LTensor_(LTensor), RTensor_(RTensor), W_(W) { // If L(R)Tensor is not present, default constructed ITensor is provided
                                size_ = 1;

                                if (LTensor_) {
                                        for (auto &I : LTensor_.inds()) {
                                                if (I.primeLevel() > 0) {
                                                        size_ *= itensor::dim(I);
                                                }
                                        }
                                }

                                if (RTensor_) {
                                        for (auto &I : RTensor_.inds()) {
                                                if (I.primeLevel() > 0) {
                                                        size_ *= itensor::dim(I);
                                                }
                                        }
                                }

                                size_ *= itensor::dim(itensor::findIndex(W_, "Site,0"));
                        }

                        size_t size() const { return size_; }

                        void product(const itensor::ITensor& phi, itensor::ITensor& phip) const {
                                phip = LTensor_ ? LTensor_ * phi : phi;
                                phip *= W_;
                                if (RTensor_) {
                                        phip *= RTensor_;
                                }
                                phip.noPrime();
                        }
        };

        class ZeroSiteH {
                private:
                        const itensor::ITensor& LTensor_, RTensor_;
                        size_t size_;
                public:
                        ZeroSiteH(const itensor::ITensor& LTensor, const itensor::ITensor& RTensor) : LTensor_(LTensor), RTensor_(RTensor) { // If L(R)Tensor is not present, default constructed ITensor is provided
                                size_ = 1;

                                if (LTensor_) {
                                        for (auto &I : LTensor_.inds()) {
                                                if (I.primeLevel() > 0) {
                                                        size_ *= itensor::dim(I);
                                                }
                                        }
                                }

                                if (RTensor_) {
                                        for (auto &I : RTensor_.inds()) {
                                                if (I.primeLevel() > 0) {
                                                        size_ *= itensor::dim(I);
                                                }
                                        }
                                }
                        }

                        size_t size() const { return size_; }

                        void product(const itensor::ITensor& phi, itensor::ITensor& phip) const {
                                phip = LTensor_ ? LTensor_ * phi : phi;
                                if (RTensor_) {
                                        phip *= RTensor_;
                                }
                                phip.noPrime();
                        }
        };
}
