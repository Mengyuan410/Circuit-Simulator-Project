#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <complex>
#include"functions.h"
#include"class.h"
#include"mode.h"





Resistor_D* D_to_R(Diode* d) {
    Resistor_D* pointer_r = new Resistor_D(d->get_anode(), d->get_cathode(), d);
    return pointer_r;
}
Resistor_M* M_to_R(Mosfet* m) {
    Resistor_M* pointer_r = new Resistor_M(m->get_drain(), m->get_source(), m);
    return pointer_r;
}
Resistor_Qbe* Q_to_Rbe(BJT* q) {
    Resistor_Qbe* pointer_r = new Resistor_Qbe(q->get_base(), q->get_emitter(), q);


    return pointer_r;
}
Resistor_Qce* Q_to_Rce(BJT* q) {

    Resistor_Qce* pointer_r = new Resistor_Qce(q->get_collector(), q->get_emitter(), q);

    return pointer_r;
}
Resistor_Qbc* Q_to_Rbc(BJT* q) {
    Resistor_Qbc* pointer_r = new Resistor_Qbc(q->get_base(), q->get_collector(), q);


    return pointer_r;
}

Vcontrolled_Isource* Q_to_Gce(BJT* q) {
    Vcontrolled_Isource* pointer_g;
    if (q->get_model() == "NPN") {
        pointer_g = new G_bjt_be(q->get_collector(), q->get_emitter(), q->get_base(), q->get_emitter(), q);
    }
    else {
        pointer_g = new G_bjt_be(q->get_emitter(), q->get_collector(), q->get_emitter(), q->get_base(), q);
    }
    return pointer_g;
}
Vcontrolled_Isource* Q_to_Gec(BJT* q) {
    Vcontrolled_Isource* pointer_g;
    if (q->get_model() == "NPN") {
        pointer_g = new G_bjt_bc(q->get_emitter(), q->get_collector(), q->get_base(), q->get_collector(), q);
    }
    else {
        pointer_g = new G_bjt_bc(q->get_collector(), q->get_emitter(), q->get_collector(), q->get_base(), q);
    }
    return pointer_g;
}
Vcontrolled_Isource* M_to_G(Mosfet* m) {
    Vcontrolled_Isource* pointer_g;
    if (m->get_model() == "NMOS") {
        pointer_g = new G_mosfet(m->get_drain(), m->get_source(), m->get_gate(), m->get_source(), m);
    }
    else {
        pointer_g = new G_mosfet(m->get_source(), m->get_drain(), m->get_gate(), m->get_source(), m);
    }
    
    return pointer_g;
}

iSource_nonlin* Q_to_Ibe(BJT* q) {
    iSource_nonlin* pointer_i;
    if (q->get_model() == "NPN") {
        pointer_i = new  I_source_Qbe(q->get_base(), q->get_emitter(), q);
    }
    else {
        pointer_i = new  I_source_Qbe(q->get_emitter(), q->get_base(), q);
    }
   
    return pointer_i;
}
iSource_nonlin* Q_to_Ibc(BJT* q) {
    iSource_nonlin* pointer_i;
    if (q->get_model() == "NPN") {
        pointer_i = new  I_source_Qbc(q->get_base(), q->get_collector(), q);
    }
    else {
        pointer_i = new  I_source_Qbc(q->get_collector(), q->get_base(), q);
    }
    return pointer_i;
}
iSource_nonlin* Q_to_Ice(BJT* q) {
    iSource_nonlin* pointer_i;
    if (q->get_model() == "NPN") {
        pointer_i = new  I_source_Qce(q->get_collector(), q->get_emitter(), q);
    }
    else {
        pointer_i = new  I_source_Qce( q->get_emitter(), q->get_collector(), q);
    }
    return pointer_i;
}
iSource_nonlin* D_to_Id(Diode* d) {
    iSource_nonlin* pointer_i = new  I_source_D(d->get_anode(), d->get_cathode());
    return pointer_i;
}
iSource_nonlin* M_to_Im(Mosfet* m) {
    iSource_nonlin* pointer_i;
    if (m->get_model() == "NMOS") {
        pointer_i = new  I_source_mosfet(m->get_drain(), m->get_source(), m);
    }
    else {
        pointer_i = new  I_source_mosfet(m->get_source(), m->get_drain(), m);
    }
    return pointer_i;
}

//resistors is all linear r. it is for dc analysis
//rlcdevices is all rlc. it is for ac analysis
//ssem_devices is rlc+all non_lin_equiv. it is for ac analysis
//large_sig_devices is r+all non_lin_equiv. it is for dc analysis
//non_lin_resis_equiv that needs to be changed for Newton Raphson method every iteration
void dc_convert(std::vector <ConductanceDevice*>& large_sig_devices, std::vector<Vcontrolled_Isource*>& G, std::vector<iSource_nonlin*>& isource_non_lin, std::vector<Resistor_non_lin_equiv*>& non_lin_resis_equiv, std::vector<ConductanceDevice*> resistors, std::vector<BJT*> BJTs, std::vector<Mosfet*> Mosfets,
    std::vector<Diode*> Diodes) {
    for (int i = 0; i < resistors.size(); i++) {
        large_sig_devices.push_back(resistors[i]);
    }

    for (int i = 0; i < Diodes.size(); i++) {
        Resistor_D* Rd = D_to_R(Diodes[i]);
        iSource_nonlin* Id = D_to_Id(Diodes[i]);
        large_sig_devices.push_back(Rd);
        non_lin_resis_equiv.push_back(Rd);
        isource_non_lin.push_back(Id);
    }
    for (int i = 0; i < BJTs.size(); i++) {
        Resistor_Qbe* Rbe = Q_to_Rbe(BJTs[i]);
        Resistor_Qce* Rce = Q_to_Rce(BJTs[i]);
        Resistor_Qbc* Rbc = Q_to_Rbc(BJTs[i]);
        Vcontrolled_Isource* Gce = Q_to_Gce(BJTs[i]);
        Vcontrolled_Isource* Gec = Q_to_Gec(BJTs[i]);
        iSource_nonlin* Ibe = Q_to_Ibe(BJTs[i]);
        iSource_nonlin* Ibc = Q_to_Ibc(BJTs[i]);
        iSource_nonlin* Ice = Q_to_Ice(BJTs[i]);
        large_sig_devices.push_back(Rbe);
        large_sig_devices.push_back(Rce);
        large_sig_devices.push_back(Rbc);
        non_lin_resis_equiv.push_back(Rbe);
        non_lin_resis_equiv.push_back(Rce);
        non_lin_resis_equiv.push_back(Rbc);
        G.push_back(Gce);
        G.push_back(Gec);
        isource_non_lin.push_back(Ibe);
        isource_non_lin.push_back(Ibc);
        isource_non_lin.push_back(Ice);
    }

    for (int i = 0; i < Mosfets.size(); i++) {
        Resistor_M* Rm = M_to_R(Mosfets[i]);
        iSource_nonlin* Im = M_to_Im(Mosfets[i]);
        Vcontrolled_Isource* Gm = M_to_G(Mosfets[i]);
        large_sig_devices.push_back(Rm);
        non_lin_resis_equiv.push_back(Rm);
        G.push_back(Gm);
        isource_non_lin.push_back(Im);
    }
}

void ac_convert(std::vector<ConductanceDevice*>& ssem_devices, const std::vector<Resistor_non_lin_equiv*> non_lin_resis_equiv, const std::vector<ConductanceDevice*> rlcdevices) {

    for (int i = 0; i < rlcdevices.size(); i++) {
        ssem_devices.push_back(rlcdevices[i]);
    }
    for (int i = 0; i < non_lin_resis_equiv.size(); i++) {
        ssem_devices.push_back(non_lin_resis_equiv[i]);
    }
}



