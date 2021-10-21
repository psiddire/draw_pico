#include "core/functions.hpp"
#include <regex>

#include "TVector2.h"

#include "core/utilities.hpp"
#include "core/config_parser.hpp"

using namespace std;

namespace Functions{

  const NamedFunc hem_veto("hem_veto",[](const Baby &b) -> NamedFunc::ScalarType{
    if(abs(b.SampleType()) == 2018) {
      if (b.SampleType()<0) {
        if (b.run() >= 319077) { 
          if(b.nel() > 0) {
            for(size_t i = 0; i < b.lep_pt()->size(); i++) {
              if (abs(b.lep_pdgid()->at(i))==13) continue;
              if(b.lep_eta()->at(i) < -1.5 && (b.lep_phi()->at(i) > -1.6 && b.lep_phi()->at(i) < -0.8)) 
                return static_cast<float>(0);
            }
          }
          for(size_t i = 0; i < b.jet_pt()->size(); i++) {
            if(Functions::IsGoodJet(b,i) && b.jet_eta()->at(i) < -1.5 && (b.jet_phi()->at(i) > -1.6 && b.jet_phi()->at(i) < -0.8)) 
              return static_cast<float>(0);
          }
        }
      } else {
        if ((b.event()%1961) < 1296) { 
          if(b.nel() > 0) {
            for(size_t i = 0; i < b.lep_pt()->size(); i++) {
              if (abs(b.lep_pdgid()->at(i))==13) continue;
              if(b.lep_eta()->at(i) < -1.5 && (b.lep_phi()->at(i) > -1.6 && b.lep_phi()->at(i) < -0.8)) 
                return static_cast<float>(0);
            }
          }
          for(size_t i = 0; i < b.jet_pt()->size(); i++) {
            if(Functions::IsGoodJet(b,i) && b.jet_eta()->at(i) < -1.5 && (b.jet_phi()->at(i) > -1.6 && b.jet_phi()->at(i) < -0.8)) 
              return static_cast<float>(0);
          }
        }
      }
    }
    return static_cast<float>(1);
  });
  
  const NamedFunc ntrub("ntrub",[](const Baby &b) -> NamedFunc::ScalarType{
    int tmp_ntrub(0);
    for (unsigned i(0); i<b.jet_pt()->size(); i++){
      if (b.jet_hflavor()->at(i)==5) tmp_ntrub++;
    }
    return tmp_ntrub;
  });

  const NamedFunc w_pileup("w_pileup",[](const Baby &b) -> NamedFunc::ScalarType{
    float weight = 1.;

    if(b.type()>0 && b.type()<1000) weight = 1; // data
    else { // mc
      vector<double> weights;
      // weights are from get_puweights.py
      if (b.SampleType()==2016) {
        weights = vector<double>({3.661e-01, 8.939e-01, 1.198e+00, 9.627e-01, 1.121e+00, 1.165e+00, 7.956e-01, 4.958e-01, 7.422e-01, 8.789e-01, 9.642e-01, 1.072e+00, 1.125e+00, 1.176e+00, 1.202e+00, 1.208e+00, 1.200e+00, 1.183e+00, 1.144e+00, 1.097e+00, 1.066e+00, 1.051e+00, 1.052e+00, 1.051e+00, 1.050e+00, 1.058e+00, 1.072e+00, 1.083e+00, 1.096e+00, 1.108e+00, 1.095e+00, 1.083e+00, 1.041e+00, 9.858e-01, 9.108e-01, 8.209e-01, 7.168e-01, 6.100e-01, 5.031e-01, 4.048e-01, 3.092e-01, 2.279e-01, 1.637e-01, 1.132e-01, 7.730e-02, 5.092e-02, 3.189e-02, 2.009e-02, 1.226e-02, 7.426e-03, 4.380e-03, 2.608e-03, 1.566e-03, 9.714e-04, 7.292e-04, 6.727e-04, 7.305e-04, 9.488e-04, 1.355e-03, 1.894e-03, 3.082e-03, 4.097e-03, 4.874e-03, 5.256e-03, 5.785e-03, 5.515e-03, 5.000e-03, 4.410e-03, 4.012e-03, 3.548e-03, 3.108e-03, 2.702e-03, 2.337e-03, 2.025e-03, 1.723e-03});
      } else if (b.SampleType()==2017) {
        weights =  vector<double>({1.835e-01, 3.933e+00, 3.471e+00, 2.492e+00, 1.625e+00, 1.515e+00, 1.288e+00, 1.278e+00, 6.158e-01, 1.452e+00, 1.498e+00, 1.487e+00, 1.331e+00, 1.164e+00, 1.079e+00, 1.054e+00, 1.081e+00, 1.129e+00, 1.166e+00, 1.189e+00, 1.213e+00, 1.238e+00, 1.260e+00, 1.271e+00, 1.273e+00, 1.271e+00, 1.271e+00, 1.267e+00, 1.274e+00, 1.252e+00, 1.221e+00, 1.170e+00, 1.109e+00, 1.038e+00, 9.694e-01, 9.119e-01, 8.668e-01, 8.345e-01, 7.884e-01, 7.508e-01, 7.592e-01, 7.932e-01, 8.583e-01, 9.588e-01, 1.094e+00, 1.256e+00, 1.420e+00, 1.496e+00, 1.531e+00, 1.462e+00, 1.337e+00, 1.154e+00, 9.506e-01, 7.497e-01, 5.698e-01, 4.107e-01, 2.900e-01, 1.987e-01, 1.376e-01, 9.664e-02, 6.923e-02, 5.086e-02, 3.844e-02, 2.996e-02, 2.409e-02, 1.712e-02, 1.248e-02, 1.077e-02, 9.597e-03, 8.812e-03, 8.310e-03, 8.018e-03, 7.885e-03, 7.875e-03, 6.337e-03, 5.334e-03, 5.444e-03, 5.584e-03, 5.744e-03, 5.916e-03, 6.090e-03, 6.257e-03, 6.408e-03, 5.095e-03, 4.228e-03, 4.259e-03, 4.266e-03, 4.247e-03, 4.203e-03, 4.137e-03, 4.052e-03, 3.950e-03, 2.943e-03, 2.296e-03, 2.200e-03, 2.101e-03, 2.002e-03, 1.902e-03, 1.804e-03});
      } else if (b.SampleType()==2018) {
        weights = vector<double>({9.687e+03, 1.327e+01, 4.377e+01, 1.869e+01, 1.251e+01, 9.042e+00, 6.579e+00, 4.873e+00, 3.629e+00, 2.763e+00, 2.224e+00, 1.900e+00, 1.702e+00, 1.579e+00, 1.504e+00, 1.465e+00, 1.450e+00, 1.451e+00, 1.459e+00, 1.464e+00, 1.458e+00, 1.437e+00, 1.401e+00, 1.353e+00, 1.301e+00, 1.250e+00, 1.205e+00, 1.167e+00, 1.137e+00, 1.115e+00, 1.100e+00, 1.089e+00, 1.082e+00, 1.077e+00, 1.073e+00, 1.069e+00, 1.063e+00, 1.052e+00, 1.036e+00, 1.012e+00, 9.794e-01, 9.373e-01, 8.862e-01, 8.272e-01, 7.618e-01, 6.923e-01, 6.211e-01, 5.505e-01, 4.827e-01, 4.194e-01, 3.618e-01, 3.105e-01, 2.659e-01, 2.276e-01, 1.952e-01, 1.679e-01, 1.450e-01, 1.259e-01, 1.097e-01, 9.590e-02, 8.401e-02, 7.361e-02, 6.442e-02, 5.623e-02, 4.889e-02, 4.231e-02, 3.643e-02, 3.122e-02, 2.662e-02, 2.262e-02, 1.915e-02, 1.618e-02, 1.363e-02, 1.147e-02, 9.633e-03, 8.075e-03, 6.750e-03, 5.623e-03, 4.661e-03, 3.838e-03, 3.134e-03, 2.531e-03, 2.016e-03, 1.578e-03, 1.208e-03, 8.991e-04, 6.465e-04, 4.458e-04, 2.929e-04, 1.824e-04, 1.094e-04, 6.258e-05, 3.428e-05, 1.807e-05, 9.221e-06, 4.578e-06, 2.221e-06, 1.058e-06, 4.964e-07, 2.308e-07});
      }
      unsigned npu_tru_mean = static_cast<unsigned> (b.npu_tru_mean());
      if (npu_tru_mean>=weights.size()) npu_tru_mean = weights.size()-1;

      weight = weights.at(npu_tru_mean);
    }

    return weight;
  });

  const NamedFunc w_syst_pileup_up("w_syst_pileup_up",[](const Baby &b) -> NamedFunc::ScalarType{
    float weight = 1.;

    if(b.type()>0 && b.type()<1000) weight = 1; // data
    else { // mc
      vector<double> weights;
      // weights are from get_puweights.py
      if (b.SampleType()==2016) {
        weights = vector<double>({3.567e-01, 7.039e-01, 1.133e+00, 8.460e-01, 1.020e+00, 1.049e+00, 7.257e-01, 3.478e-01, 5.006e-01, 6.030e-01, 6.322e-01, 7.328e-01, 8.278e-01, 9.125e-01, 9.599e-01, 9.885e-01, 1.024e+00, 1.053e+00, 1.051e+00, 1.027e+00, 1.006e+00, 9.980e-01, 1.015e+00, 1.038e+00, 1.058e+00, 1.085e+00, 1.121e+00, 1.155e+00, 1.193e+00, 1.231e+00, 1.246e+00, 1.268e+00, 1.259e+00, 1.234e+00, 1.182e+00, 1.105e+00, 1.002e+00, 8.894e-01, 7.691e-01, 6.539e-01, 5.325e-01, 4.226e-01, 3.298e-01, 2.499e-01, 1.886e-01, 1.382e-01, 9.684e-02, 6.865e-02, 4.735e-02, 3.249e-02, 2.171e-02, 1.453e-02, 9.585e-03, 6.168e-03, 4.278e-03, 3.071e-03, 2.242e-03, 1.923e-03, 2.023e-03, 2.397e-03, 3.652e-03, 4.779e-03, 5.725e-03, 6.274e-03, 7.045e-03, 6.866e-03, 6.372e-03, 5.755e-03, 5.366e-03, 4.866e-03, 4.374e-03, 3.905e-03, 3.469e-03, 3.090e-03, 2.704e-03});
      } else if (b.SampleType()==2017) {
        weights = vector<double>({1.777e-01, 3.395e+00, 2.840e+00, 2.506e+00, 1.399e+00, 1.342e+00, 1.244e+00, 1.169e+00, 4.957e-01, 9.402e-01, 1.066e+00, 1.034e+00, 1.009e+00, 9.005e-01, 8.268e-01, 7.968e-01, 8.116e-01, 8.686e-01, 9.339e-01, 9.858e-01, 1.034e+00, 1.085e+00, 1.128e+00, 1.155e+00, 1.165e+00, 1.167e+00, 1.172e+00, 1.182e+00, 1.210e+00, 1.212e+00, 1.207e+00, 1.179e+00, 1.140e+00, 1.089e+00, 1.039e+00, 9.943e-01, 9.571e-01, 9.283e-01, 8.793e-01, 8.335e-01, 8.292e-01, 8.394e-01, 8.679e-01, 9.197e-01, 1.000e+00, 1.112e+00, 1.250e+00, 1.347e+00, 1.448e+00, 1.488e+00, 1.490e+00, 1.429e+00, 1.316e+00, 1.165e+00, 9.945e-01, 8.024e-01, 6.303e-01, 4.766e-01, 3.612e-01, 2.748e-01, 2.113e-01, 1.650e-01, 1.315e-01, 1.070e-01, 8.913e-02, 6.499e-02, 4.817e-02, 4.188e-02, 3.730e-02, 3.399e-02, 3.168e-02, 3.016e-02, 2.930e-02, 2.901e-02, 2.327e-02, 1.963e-02, 2.021e-02, 2.104e-02, 2.208e-02, 2.331e-02, 2.469e-02, 2.620e-02, 2.781e-02, 2.297e-02, 1.986e-02, 2.090e-02, 2.191e-02, 2.288e-02, 2.381e-02, 2.470e-02, 2.553e-02, 2.632e-02, 2.077e-02, 1.719e-02, 1.750e-02, 1.779e-02, 1.806e-02, 1.832e-02, 1.856e-02});
      } else if (b.SampleType()==2018) {
        weights = vector<double>({9.001e+03, 1.158e+01, 3.735e+01, 1.605e+01, 1.078e+01, 7.826e+00, 5.625e+00, 4.121e+00, 3.054e+00, 2.298e+00, 1.815e+00, 1.527e+00, 1.356e+00, 1.255e+00, 1.196e+00, 1.167e+00, 1.160e+00, 1.168e+00, 1.187e+00, 1.208e+00, 1.226e+00, 1.234e+00, 1.230e+00, 1.211e+00, 1.182e+00, 1.147e+00, 1.112e+00, 1.080e+00, 1.053e+00, 1.033e+00, 1.020e+00, 1.013e+00, 1.010e+00, 1.012e+00, 1.018e+00, 1.026e+00, 1.035e+00, 1.044e+00, 1.052e+00, 1.056e+00, 1.054e+00, 1.046e+00, 1.029e+00, 1.003e+00, 9.684e-01, 9.245e-01, 8.730e-01, 8.154e-01, 7.539e-01, 6.903e-01, 6.268e-01, 5.652e-01, 5.069e-01, 4.530e-01, 4.041e-01, 3.605e-01, 3.220e-01, 2.884e-01, 2.593e-01, 2.340e-01, 2.120e-01, 1.926e-01, 1.754e-01, 1.598e-01, 1.456e-01, 1.324e-01, 1.202e-01, 1.087e-01, 9.800e-02, 8.804e-02, 7.884e-02, 7.040e-02, 6.270e-02, 5.572e-02, 4.943e-02, 4.376e-02, 3.867e-02, 3.409e-02, 2.995e-02, 2.621e-02, 2.280e-02, 1.968e-02, 1.681e-02, 1.415e-02, 1.169e-02, 9.426e-03, 7.364e-03, 5.535e-03, 3.975e-03, 2.714e-03, 1.790e-03, 1.129e-03, 6.830e-04, 3.989e-04, 2.259e-04, 1.248e-04, 6.746e-05, 3.584e-05, 1.877e-05, 9.714e-06});
      }
      unsigned npu_tru_mean = static_cast<unsigned> (b.npu_tru_mean());
      if (npu_tru_mean>=weights.size()) npu_tru_mean = weights.size()-1;

      weight = weights.at(npu_tru_mean);
    }

    return weight;
  });

  const NamedFunc w_syst_pileup_down("w_syst_pileup_down",[](const Baby &b) -> NamedFunc::ScalarType{
    float weight = 1.;

    if(b.type()>0 && b.type()<1000) weight = 1; // data
    else { // mc
      vector<double> weights;
      // weights are from get_puweights.py
      if (b.SampleType()==2016) {
        weights = vector<double>({3.793e-01, 1.141e+00, 1.260e+00, 1.099e+00, 1.250e+00, 1.281e+00, 9.202e-01, 7.677e-01, 1.093e+00, 1.337e+00, 1.486e+00, 1.528e+00, 1.498e+00, 1.501e+00, 1.497e+00, 1.444e+00, 1.368e+00, 1.299e+00, 1.227e+00, 1.166e+00, 1.125e+00, 1.091e+00, 1.064e+00, 1.040e+00, 1.019e+00, 1.006e+00, 9.970e-01, 9.849e-01, 9.728e-01, 9.565e-01, 9.147e-01, 8.722e-01, 8.071e-01, 7.343e-01, 6.510e-01, 5.614e-01, 4.664e-01, 3.745e-01, 2.885e-01, 2.145e-01, 1.497e-01, 9.989e-02, 6.434e-02, 3.959e-02, 2.389e-02, 1.382e-02, 7.566e-03, 4.156e-03, 2.215e-03, 1.187e-03, 6.427e-04, 3.841e-04, 2.728e-04, 2.434e-04, 2.896e-04, 3.958e-04, 5.440e-04, 7.828e-04, 1.151e-03, 1.604e-03, 2.568e-03, 3.340e-03, 3.880e-03, 4.079e-03, 4.373e-03, 4.058e-03, 3.579e-03, 3.068e-03, 2.712e-03, 2.327e-03, 1.978e-03, 1.667e-03, 1.397e-03, 1.172e-03, 9.649e-04});
      } else if (b.SampleType()==2017) {
        weights = vector<double>({1.919e-01, 4.353e+00, 4.534e+00, 2.409e+00, 1.874e+00, 1.690e+00, 1.344e+00, 1.433e+00, 8.913e-01, 2.162e+00, 2.191e+00, 2.070e+00, 1.742e+00, 1.533e+00, 1.440e+00, 1.424e+00, 1.433e+00, 1.434e+00, 1.425e+00, 1.408e+00, 1.394e+00, 1.390e+00, 1.391e+00, 1.393e+00, 1.392e+00, 1.381e+00, 1.361e+00, 1.329e+00, 1.308e+00, 1.257e+00, 1.201e+00, 1.125e+00, 1.042e+00, 9.550e-01, 8.776e-01, 8.173e-01, 7.734e-01, 7.466e-01, 7.160e-01, 7.041e-01, 7.473e-01, 8.269e-01, 9.431e-01, 1.089e+00, 1.248e+00, 1.392e+00, 1.483e+00, 1.436e+00, 1.323e+00, 1.124e+00, 9.070e-01, 6.896e-01, 5.007e-01, 3.501e-01, 2.378e-01, 1.547e-01, 9.969e-02, 6.300e-02, 4.070e-02, 2.693e-02, 1.836e-02, 1.297e-02, 9.519e-03, 7.283e-03, 5.810e-03, 4.134e-03, 3.040e-03, 2.659e-03, 2.404e-03, 2.235e-03, 2.124e-03, 2.052e-03, 2.007e-03, 1.979e-03, 1.563e-03, 1.283e-03, 1.271e-03, 1.260e-03, 1.248e-03, 1.234e-03, 1.215e-03, 1.191e-03, 1.160e-03, 8.750e-04, 6.873e-04, 6.539e-04, 6.171e-04, 5.777e-04, 5.367e-04, 4.949e-04, 4.532e-04, 4.125e-04, 2.864e-04, 2.079e-04, 1.851e-04, 1.641e-04, 1.448e-04, 1.273e-04, 1.116e-04});
      } else if (b.SampleType()==2018) {
        weights = vector<double>({1.048e+04, 1.536e+01, 5.166e+01, 2.189e+01, 1.460e+01, 1.054e+01, 7.785e+00, 5.822e+00, 4.363e+00, 3.381e+00, 2.776e+00, 2.402e+00, 2.164e+00, 2.008e+00, 1.910e+00, 1.853e+00, 1.822e+00, 1.805e+00, 1.787e+00, 1.757e+00, 1.711e+00, 1.647e+00, 1.572e+00, 1.494e+00, 1.422e+00, 1.358e+00, 1.306e+00, 1.265e+00, 1.232e+00, 1.206e+00, 1.184e+00, 1.166e+00, 1.148e+00, 1.131e+00, 1.110e+00, 1.085e+00, 1.054e+00, 1.014e+00, 9.650e-01, 9.068e-01, 8.402e-01, 7.668e-01, 6.889e-01, 6.091e-01, 5.303e-01, 4.549e-01, 3.851e-01, 3.223e-01, 2.673e-01, 2.202e-01, 1.808e-01, 1.483e-01, 1.218e-01, 1.003e-01, 8.303e-02, 6.903e-02, 5.761e-02, 4.823e-02, 4.042e-02, 3.386e-02, 2.831e-02, 2.358e-02, 1.954e-02, 1.610e-02, 1.319e-02, 1.074e-02, 8.700e-03, 7.011e-03, 5.627e-03, 4.500e-03, 3.588e-03, 2.854e-03, 2.265e-03, 1.792e-03, 1.413e-03, 1.110e-03, 8.663e-04, 6.717e-04, 5.164e-04, 3.929e-04, 2.953e-04, 2.186e-04, 1.591e-04, 1.133e-04, 7.865e-05, 5.292e-05, 3.429e-05, 2.125e-05, 1.250e-05, 6.957e-06, 3.719e-06, 1.892e-06, 9.197e-07, 4.300e-07, 1.947e-07, 8.630e-08, 3.802e-08, 1.715e-08, 8.357e-09, 4.754e-09});
      }
      unsigned npu_tru_mean = static_cast<unsigned> (b.npu_tru_mean());
      if (npu_tru_mean>=weights.size()) npu_tru_mean = weights.size()-1;

      weight = weights.at(npu_tru_mean);
    }

    return weight;
  });

  const NamedFunc SampleType("SampleType",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.SampleType();
  });


  bool IsGoodJet(const Baby &b, size_t ijet){
    return ijet<b.jet_pt()->size()
      && b.jet_pt()->at(ijet) > 30.
      && fabs(b.jet_eta()->at(ijet))<2.4
      && !b.jet_islep()->at(ijet);
  }

  NamedFunc MismeasurementCorrection(const string &file_path,
                                     const string &mismeas_scenario,
                                     Variation variation){
    ConfigParser cp;
    cp.Load(file_path, mismeas_scenario);

    string name = "w_mismeas";
    NamedFunc reweight_cut = cp.GetOpt<std::string>("reweight_cut");
    double wgt = 1.;
    switch(variation){
    case Variation::central:
      wgt = cp.GetOpt<double>("w_central");
      break;
    case Variation::up:
      wgt = cp.GetOpt<double>("w_up");
      break;
    case Variation::down:
      wgt = cp.GetOpt<double>("w_down");
      break;
    default:
      wgt = 1.;
      ERROR("Invalid variation: "+to_string(static_cast<int>(variation)));
      break;
    }

    if(reweight_cut.IsScalar()){
      return NamedFunc(name, [reweight_cut, wgt](const Baby &b){
          return reweight_cut.GetScalar(b) ? wgt : 1.;
        });
    }else{
      return NamedFunc(name, [reweight_cut, wgt](const Baby &b){
          NamedFunc::VectorType cuts = reweight_cut.GetVector(b);
          NamedFunc::VectorType wgts(cuts.size());
          for(size_t i = 0; i < cuts.size(); ++i){
            wgts.at(i) = cuts.at(i) ? wgt : 1.;
          }
          return wgts;
        });
    }
  }

  NamedFunc MismeasurementWeight(const std::string &file_path,
                                 const std::string &mismeas_scenario){
    ConfigParser cp;
    cp.Load(file_path, mismeas_scenario);
    NamedFunc cut = cp.GetOpt<std::string>("mismeas_cut");
    NamedFunc wgt = cp.GetOpt<std::string>("mismeas_wgt");
    string name = "w_sys_mm";
    if(cut.IsScalar()){
      if(wgt.IsScalar()){
        return NamedFunc(name, [cut, wgt](const Baby &b){
            return cut.GetScalar(b) ? wgt.GetScalar(b) : 1.;
          });
      }else{
        return NamedFunc(name, [cut, wgt](const Baby &b){
            NamedFunc::VectorType wgts = wgt.GetVector(b);
            return cut.GetScalar(b) ? wgts : NamedFunc::VectorType(wgts.size(), 1.);
          });
      }
    }else{
      if(wgt.IsScalar()){
        return NamedFunc(name, [cut, wgt](const Baby &b){
            NamedFunc::VectorType cuts = cut.GetVector(b);
            NamedFunc::VectorType wgts(cuts.size());
            NamedFunc::ScalarType scalar_wgt = wgt.GetScalar(b);
            for(size_t i = 0; i < wgts.size(); ++i){
              wgts.at(i) = cuts.at(i) ? scalar_wgt : 1.;
            }
            return wgts;
          });
      }else{
        return NamedFunc(name, [cut, wgt](const Baby &b){
            NamedFunc::VectorType cuts = cut.GetVector(b);
            NamedFunc::VectorType wgts = wgt.GetVector(b);
            size_t min_size = min(cuts.size(), wgts.size());
            NamedFunc::VectorType out(min_size);
            for(size_t i = 0; i < min_size; ++i){
              out.at(i) = cuts.at(i) ? wgts.at(i) : 1.;
            }
            return out;
          });
      }
    }
  }

  string regSearch(string inString, string const & regExPattern) {
      std::smatch match;
      std::regex pattern(regExPattern);
      std::regex_search (inString, match, pattern);
      return match.str();
  }
  pair<int, int> getSignalMassValues(string filename) {
    pair<int, int> nlsp_lsp_mass;
    string baseFilename = filename.substr(filename.find_last_of("/\\")+1);
    string datasetName = regSearch(baseFilename, "[A-Z].*|ttHTobb.*");
    if (datasetName.find("SMS-TChiHH") == string::npos) return {-1,-1};
    string NLSPMassString = regSearch(datasetName, "mChi-.*?_").substr(5);
    string LSPMassString = regSearch(datasetName, "mLSP-.*?_").substr(5);
    //cout<<NLSPMassString<<" "<<LSPMassString<<endl;
    return {stoi(NLSPMassString),stoi(LSPMassString)};
  }
  bool findStringIC(const std::string & strHaystack, const std::string & strNeedle)
  {
    auto it = std::search(
      strHaystack.begin(), strHaystack.end(),
      strNeedle.begin(),   strNeedle.end(),
      [](char ch1, char ch2) { return std::toupper(ch1) == std::toupper(ch2); }
    );
    return (it != strHaystack.end() );
  }

  bool regionCut(const Baby & b, int regionIndex) {
    bool inRegion = false;
    vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > > * eventNumberData = static_cast<vector<map<Long64_t, set< tuple<string, int, int, int, int, int, int> > > > *> (b.EventVetoData());
    // 0: sampleType, 1: year, 2: run, 3: lumiblock, 4: met, 5: nlsp_mass, 6: lsp_mass
    // Check event nuber
    if ((*eventNumberData)[regionIndex].count(b.event())==1) {
      for (auto const & it : (*eventNumberData)[regionIndex][b.event()]) {
        // Check year
        if (get<1>(it) != abs(b.SampleType())) continue;
        // Check run
        if (get<2>(it) != b.run()) continue;
        bool isSignal = ((*b.FileNames().begin()).find("TChiHH") != string::npos) ? true : false;
        // Check met out of 20%
        if (isSignal) {if (get<4>(it) < b.met()*0.8 || get<4>(it) > b.met()*1.2) continue;}
        else if (get<4>(it) != round(b.met())) continue;
        pair<int, int> signal_mass = getSignalMassValues(*b.FileNames().begin());
        // Check nlsp_mass
        if (signal_mass.first != get<5>(it)) continue;
        // Check lsp_mass
        int regionLsp = get<6>(it)==1? 0:get<6>(it);
        if (signal_mass.second != regionLsp) continue;
        // Check sample name
        //if (get<5>(it) != -1) cout<<*b.FileNames().begin()<<" "<<get<0>(it)<<endl;
        if ((*b.FileNames().begin()).find("TChiHH") != string::npos && get<0>(it).find("TChiHH") != string::npos)  {
          // Need to reject if 1D is compared with 2D
          int babyDimension = (*b.FileNames().begin()).find("HToBB_2D") == string::npos ? 1 : 2;
          int eventListDimension = get<0>(it).find("2D") == string::npos ? 1 : 2;
          if (babyDimension != eventListDimension) continue;
          //if ((*b.FileNames().begin()).find("HToBB_2D") != get<0>(it).find("2D")) continue;
        } else if (!findStringIC(*b.FileNames().begin(), get<0>(it))) continue;
        //if (isSignal) cout<<"event: "<<b.event()<<" => sampleType: "<<get<0>(it)<<" year: "<<get<1>(it)<<" run: "<<get<2>(it)<<" lumiblock: "<<get<3>(it)<<" met: "<<get<4>(it)<<" nlsp_mass: "<<get<5>(it)<<" lsp_mass: "<<get<6>(it)<<endl;
  
        //if (b.met()<300) cout<<get<0>(it)<<" "<<get<1>(it)<<" "<<get<2>(it)<<" "<<get<3>(it)<<" "<<b.event()<<endl;
        //cout<<"event: "<<b.event()<<" => sampleType: "<<get<0>(it)<<" year: "<<get<1>(it)<<" run: "<<get<2>(it)<<" lumiblock: "<<get<3>(it)<<" met: "<<get<4>(it)<<" nlsp_mass: "<<get<5>(it)<<" lsp_mass: "<<get<6>(it)<<endl;
  
        inRegion = true;
        break;
      }
    }
    return inRegion;
  }
  const NamedFunc boostSignalRegion("boostSignalRegion",[](const Baby &b) -> NamedFunc::ScalarType{
    return regionCut(b, 0);
  });
  const NamedFunc boostControlRegion("boostControlRegion",[](const Baby &b) -> NamedFunc::ScalarType{
    return regionCut(b, 1);
  });

}
