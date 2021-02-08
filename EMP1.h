//
// Created by dav0 on 2/1/21.
//

#ifndef UNTITLED_EMP1_H
#define UNTITLED_EMP1_H

#include <vector>
#include <random>

bool read_csv_file(const char *csv_name, std::vector<std::vector<std::string>> &vec_str);

class EMP1 {
public:
    EMP1(int N, const char *pv_csv_file_name, const char *i0_csv_file_name);

    ~EMP1();

    std::tuple<int, int, int, double, double> EMP1_DRA(double pnk, double sk);

    void Simulate(double s0);

    /**
   * @brief Generate PV Power values
   *
   * Computes all the required quantities.
   *
   */
    void GeneratePv();

    /**
    * @brief Same as GeneratePV() but loads the C_PV parameter
    *
    *
    *
    */
    void GeneratePv(double new_cpv);

    double getCpv() const;
    void setCpv(double cpv);
    double getCb() const;
    void setCb(double cb);
    double getSTilde() const;
    void setSTilde(double sTilde);

private:
    EMP1(const EMP1 &);

    void operator=(const EMP1 &);

    /**
   * @brief Load the data from the CSV files
   *
   * Read the data from the plain comma separated text files.
   *
   */
    void LoadData();

    void GenerateKtSamples();

    void GenerateIrradation();

    void GenerateLoad();


    std::vector<double> kt_;
    std::vector<double> i_;
    std::vector<double> i0_data_;
    std::vector<double> pl_data_;
    std::vector<double> rk_;
    std::vector<double> pL_;
    std::vector<double> pPv_;
    const char *pv_csv_;
    const char *i0_csv_;

    std::vector<std::vector<std::string>> data_string;
    std::vector<std::tuple<int, int, int>> m_d_h_; //month day hour.
    std::vector<std::tuple<int, int, int, double, double>> sigma_u_vec_;

    /** Params Horizon Length */
    int N_;
    /** Params: Reflectance of the panel */
    double reflectance_{0.6};
    /** Params latitud */
    double latitude_{40.};
    /** Params tilt angle */
    double tilt_{32.};
    /** Params reflecting rho */
    double rating_f_{0.95};

    double CPV_;
    double CB_;
    double s_tilde_;

    double p_d_;
    double s_lb_;
    double s_ub_;

    std::normal_distribution<double> *nd_gen_f_;
};

double id_i_from_kt(double kt);

#endif //UNTITLED_EMP1_H
