//
// Created by dav0 on 2/1/21.
//

#include <cmath>

#include <chrono>
#include <iostream>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <tuple>
#include <iomanip>
#include <map>
#include "EMP1.h"

#define PI 3.14159265


EMP1::EMP1(int N, const char *pv_csv_file_name, const char *i0_csv_file_name)
        : N_(N),
          reflectance_(0.6),
          pv_csv_(pv_csv_file_name),
          i0_csv_(i0_csv_file_name),
          CPV_(20.),
          CB_(350.),
          s_tilde_(0.9),
          p_d_(40.),
          s_lb_(0.3),
          s_ub_(0.99)
          {

    nd_gen_f_ = new std::normal_distribution<double>(0., 0.003);
    LoadData();
    GenerateKtSamples();
    GenerateIrradation();
    GenerateLoad();
    GeneratePv();

}

EMP1::~EMP1() = default;

void EMP1::GenerateKtSamples() {
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine dre(seed);
    double k = 0.5;  // initial k (clearness index)
    for (int i = 0; i < N_; i++) {
        double e = (*nd_gen_f_)(dre);
        k = k * 1. + e;
        k = k > 0. ? k : 0.;   // enforce k \in [0, 1]
        k = (k <= 1.) ? k : 1.;
        kt_.push_back(k);
    }
    //
}

void EMP1::GenerateIrradation() {
    double id_i;
    double ib_i;
    double rb_ = 1.;
    double I_T;
    double hour_angle;
    double declination_angle;
    int day_of_year;
    std::map<int, int> days_accumulative;
    days_accumulative.insert(std::pair<int, int>(1, 0));
    days_accumulative.insert(std::pair<int, int>(2, 31));
    days_accumulative.insert(std::pair<int, int>(3, days_accumulative[2] + 28));
    days_accumulative.insert(std::pair<int, int>(4, days_accumulative[3] + 31));
    days_accumulative.insert(std::pair<int, int>(5, days_accumulative[4] + 30));
    days_accumulative.insert(std::pair<int, int>(6, days_accumulative[5] + 31));
    days_accumulative.insert(std::pair<int, int>(7, days_accumulative[6] + 30));
    days_accumulative.insert(std::pair<int, int>(8, days_accumulative[7] + 31));
    days_accumulative.insert(std::pair<int, int>(9, days_accumulative[8] + 31));
    days_accumulative.insert(std::pair<int, int>(10, days_accumulative[9] + 30));
    days_accumulative.insert(std::pair<int, int>(11, days_accumulative[10] + 31));
    days_accumulative.insert(std::pair<int, int>(12, days_accumulative[11] + 30));


    int j = 0;
    for (const auto &i : kt_) {
        id_i = id_i_from_kt(i);
        ib_i = 1. - id_i;
        std::cout << std::setw(3) << std::get<1>(m_d_h_[j]) << " " << std::setw(4) << std::get<2>(m_d_h_[j]);


        day_of_year = days_accumulative[std::get<0>(m_d_h_[j])] + std::get<1>(m_d_h_[j]);
        std::cout << "\t" << std::setw(4) << day_of_year << "\t";
        std::cout << std::setw(8) << i << "\t";
        std::cout << std::setw(8) << id_i << "\t" << std::setw(6) << ib_i << "\t";

        hour_angle = (std::get<2>(m_d_h_[j]) - 12) * 15.;
        hour_angle = static_cast<double>(hour_angle);
        declination_angle = 23.45 * sin(360. * (284. + day_of_year) / 365.);
        if (declination_angle < -23.45 || declination_angle > 23.45) {
            std::cout << "\nError: declination angle out of bounds.\n";
        }

        // Ratio of beam radiation on tilted surface to that on horizontal surface:
        rb_ = (cos((latitude_ - tilt_) * PI / 180.) * cos(declination_angle * PI / 180.) * cos(hour_angle * PI / 180.) +
               sin((latitude_ - tilt_) * PI / 180.) * sin(declination_angle * PI / 180.)) /
              (cos(latitude_ * PI / 180.) * cos(declination_angle * PI / 180.) * cos(hour_angle * PI / 180.) +
               sin(latitude_ * PI / 180.) * sin(declination_angle * PI / 180.));

        std::cout << std::setw(6) << hour_angle << "\t";
        std::cout << std::setw(6) << cos(hour_angle * PI / 180.) << "\t";
        std::cout << std::setw(9) << declination_angle << "\t";
        std::cout << std::setw(8) << rb_ << "\t";

        I_T = i0_data_[j] * (ib_i * rb_ + id_i * (1. + cos(tilt_)) / 2. + reflectance_ * (1. - cos(tilt_)) / 2.);
        I_T = I_T < 0.0 ? 0.0 : I_T;
        i_.push_back(I_T);

        std::cout << std::setw(8) << i0_data_[j] << "\t";
        std::cout << std::setw(8) << I_T << "\n";

        j++;
    }
}


void EMP1::LoadData() {
    bool stat = true;
    char *end_ptr;
    stat = read_csv_file(pv_csv_, data_string);
    if (stat) {
        throw "Failure reading pv csv file.\n";
    }
    for (const auto &i : data_string) {
        std::tuple<int, int, int> m_d_h;
        std::get<0>(m_d_h) = stoi(i[0]);
        std::get<1>(m_d_h) = stoi(i[1]);
        std::get<2>(m_d_h) = stoi(i[2]);
        m_d_h_.push_back(m_d_h);
        double d = std::strtod(i[3].c_str(), &end_ptr);
        pl_data_.push_back(d);
    }
    data_string.clear();
    stat = read_csv_file(i0_csv_, data_string);
    if (stat) {
        throw "Failure reading i0 csv file.\n";
    }

    for (const auto &i: data_string) {
        double d = std::strtod(i.back().c_str(), &end_ptr);
        i0_data_.push_back(d);
    }
    /*
    for(int i = 0; i < i0_data_.size(); i++){
        std::cout << i << "\t";
        std::cout << std::setw(5) << pl_data_[i];
        std::cout << "\t";
        std::cout << std::setw(5) << i0_data_[i];
        std::cout << "\n";
    }
    */
    std::cout << "Data info:" << "\n";
    std::cout << "PV size = " << pl_data_.size() << "\n";
    std::cout << "I0 size = " << i0_data_.size() << "\n";
    if (!(pl_data_.size() == i0_data_.size())) {
        std::cout << "The sizes of the data vectors is inconsistent." << "\n";
    }
    N_ = pl_data_.size() > N_ ? pl_data_.size() : static_cast<unsigned long>(N_);

}

double EMP1::getCpv() const {
    return CPV_;
}

void EMP1::setCpv(double cpv) {
    CPV_ = cpv;
}

double EMP1::getCb() const {
    return CB_;
}

void EMP1::setCb(double cb) {
    CB_ = cb;
}

double EMP1::getSTilde() const {
    return s_tilde_;
}

void EMP1::setSTilde(double sTilde) {
    s_tilde_ = sTilde;
}

void EMP1::GenerateLoad() {
    std::normal_distribution<double> nrm_dis_load(0.0, 1);
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine  dre(seed);
    for (int i = 0; i < N_; i++){
        double r = nrm_dis_load(dre);
        r = r < (1E+03) ? r : 1E+03;
        r = r > (-1E+03) ? r : -1E+03;
        rk_.push_back(r); // for debugging purposes.
        pL_.push_back(pl_data_[i] + r);

        //std::cout << r << "\t" << pl_data_[i] << "\n";
    }


}

void EMP1::GeneratePv() {
    double pv;
    for(const auto& i : i_){
        pv = rating_f_ * CPV_ * ((i / 1000.) / 1.);
        pPv_.push_back(pv);
    }
}

void EMP1::GeneratePv(double new_cpv) {
    setCpv(new_cpv);
    GeneratePv();
}

std::tuple<int, int, int, double, double> EMP1::EMP1_DRA(double pnk, double sk) {
    double s1, s2, s_next, u_k;
    int sig_1, sig_2, sig_3, dk;

    s1 = sk + pnk/CB_;
    sig_1 = (s1 <= s_tilde_) ? 1 : -1;
    dk = s1 <= s_tilde_ ? 1 : 0;
    s2 = sk + (dk * p_d_ + pnk)/CB_;
    sig_2 = s2 <= s_lb_ ? 1: -1;
    if (sig_2 == -1){
        sig_3 = s2 <= s_ub_ ? 1 : -1;
        s_next = s2 <= s_ub_ ? s2: s_ub_;
        u_k = 0;
    } else {
        if (s2 <= s_ub_){
            sig_3 = 1;
            s_next = s_lb_;
            u_k = (s_lb_ - s2) * CB_;
        } else {
            throw "whoops";
        }
    }
    return std::make_tuple(sig_1, sig_2, sig_3, u_k, s_next);

}

void EMP1::Simulate(double s0) {
    int i;
    double pn;
    // state-
    //of-charge (SOC) of the battery
    double sk = s0;
    size_t hW[7] = {
            std::string("Hour").size(),
            std::string("pn Value").size() + 6,
            std::string("s Value").size() + 4,
            std::string("Sigma1").size(),
            std::string("Sigma2").size(),
            std::string("Sigma3").size(),
            std::string("u Value").size() + 4
    };
    std::cout << std::left << "Hour" << " ";
    std::cout << "pn Value      " << " ";
    std::cout << "s Value    " << " ";
    std::cout << std::right << "Sigma1" << " ";
    std::cout << "Sigma2" << " ";
    std::cout << "Sigma3" << " ";
    std::cout << "u Value" << "\n";
    for (i = 0; i < N_; i++){
        pn = pPv_[i] - pL_[i];
        std::tuple<int, int, int, double, double> l_res = EMP1_DRA(pn, sk);
        sigma_u_vec_.push_back(l_res);
        sk = std::get<4>(l_res);
        std::cout << std::left << std::setw(hW[0]) << i << " ";
        std::cout << std::setw(hW[1]) << pn << " ";
        std::cout << std::setw(hW[2]) << sk << " ";
        std::cout << std::right << std::setw(hW[3]) << std::get<0>(l_res) << " ";
        std::cout << std::setw(hW[4]) << std::get<1>(l_res) << " ";
        std::cout << std::setw(hW[5]) << std::get<2>(l_res) << " ";
        std::cout << std::setw(hW[6]) << std::get<3>(l_res) << "\n";

    }

}

double id_i_from_kt(double kt) {
    double id_i;
    if (kt <= 0.22) {
        id_i = 1. - 0.09 * kt;
    } else if (0.22 < kt <= 0.8) {
        id_i = 0.9811 - 0.1604 * kt + pow(kt, 2.);
    } else {
        id_i = 0.165;
    }
    return id_i;
}


bool read_csv_file(const char *csv_name, std::vector<std::vector<std::string>> &vec_str) {
    std::fstream csv_file;
    csv_file.open(csv_name, std::ios::in);
    /*
     * If not open, return true.
     */
    if (!csv_file.is_open()) {
        return true;
    }

    std::string line;
    std::vector<std::vector<std::string>> data;


    if (csv_file.is_open()) {
        while (getline(csv_file, line)) {
            std::vector<std::string> string_v;
            boost::tokenizer<boost::escaped_list_separator<char>> token(line);
            string_v.assign(token.begin(), token.end());
            data.push_back(string_v);
        }
    }
    csv_file.close();
    vec_str = data;
    return false;
}