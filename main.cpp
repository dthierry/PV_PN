#include <iostream>
#include <fstream>
#include <string>
#include <boost/tokenizer.hpp>
#include <algorithm>
#include <iterator>
#include "EMP1.h"


int main() {

    /*
    for (const auto& i : mydata){
        copy(i.begin(), i.end(), std::ostream_iterator<std::string>(std::cout, "|"));
        std::cout << "\n";
    }
     */

    /*
     * Begin
     */

    /** Create emp1 object */
    EMP1 emp1(10, "e_loads_Texas_MDHPD_kWh.csv", "irradiance_MDH_G_I_Wm2.csv");
    /** Generate PhotoVoltaic Power */
    emp1.GeneratePv(20.);
    /** Set design capacity */
    emp1.setCb(350.);
    /** As a function of initial s */
    emp1.Simulate(0.6);

    return 0;
}
