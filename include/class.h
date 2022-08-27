#ifndef classes
#define classes

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <complex>
#include <iterator>
#include <cmath>
#include "mode.h"

static double Vt = 0.025;
static double Is = std::pow(10, -14);


class Node {
public:
    Node(int n) : node_orig(n) {
        node_dc = n;
        node_ac = n;
        voltage_dc = 0;
        voltage_ac = 0;
    }
    int get_node_orig() const { return node_orig; }
    int get_node_dc() const { return node_dc; }
    int get_node_ac() const { return node_ac; }
    std::complex<double> get_voltage_dc() const { return voltage_dc; }
    std::complex<double> get_voltage_ac() const { return voltage_ac; }
    void change_node_dc(int new_node) { node_dc = new_node; }
    void change_node_ac(int new_node) { node_ac = new_node; }
    void change_voltage_dc(std::complex<double> voltage) { voltage_dc = voltage; }
    void change_voltage_ac(std::complex<double> voltage) { voltage_ac = voltage; }
private:
    int node_orig;
    int node_dc;
    int node_ac;
    std::complex<double> voltage_dc;
    std::complex<double> voltage_ac;
};

class Diode {
public:
    // constructor with initializer list
    // current flows from anode to cathod
    Diode(Node* n1, Node* n2, std::string m) : anode(n1), cathode(n2), model(m) { }

    Node* get_anode() const {
        return anode;
    }

    Node* get_cathode() const {
        return cathode;
    }

private:
    Node* anode;
    Node* cathode;
    std::string model;
};

class BJT {
public:
    // constructor with initializer list
    BJT(Node* c, Node* b, Node* e, std::string m) : collector(c), base(b), emitter(e), model(m) { }

    Node* get_collector() const {
        return collector;
    }

    Node* get_base() const {
        return base;
    }

    Node* get_emitter() const {
        return emitter;
    }

    std::string get_model() const {
        return model;
    }

    double get_beta_forward() const {
        // npn: 200  pnp: 250
        if (model == "NPN") {
            return 200;
        }
        else if (model == "PNP") {
            return 250;
        }
        else {
            return false;
        }
    }

    double get_beta_reverse() const {
        // npn: 3  pnp: 3
        return 3;
    }

    double get_VA() const {
        // npn: 100  pnp: 120
        if (model == "NPN") {
            return 100;
        }
        else if (model == "PNP") {
            return 120;
        }
        else {
            return false;
        }
    }

private:
    Node* collector;
    Node* base;
    Node* emitter;
    std::string model;
};

class Mosfet {
public:
    // constructor with initializer list
    // m is "NMOS" or "PMOS"
    Mosfet(Node* d, Node* g, Node* s, std::string m) : drain(d), gate(g), source(s), model(m) { }

    Node* get_drain() const {
        return drain;
    }

    Node* get_gate() const {
        return gate;
    }

    Node* get_source() const {
        return source;
    }

    double get_k() const {
        return 0.0025;
    }

    double get_VA() const {
        return 100;
    }

    std::string get_model() const {
        return model;
    }

    double get_Vt_mosfet() const {
        if (model == "NMOS") {
            return 2;
        }
        else if (model == "PMOS") {
            return -2;
        }
        else {
            return false;
        }
    }

private:
    Node* drain;
    Node* gate;
    Node* source;
    std::string model;
};


// Parent class: ConductanceDevice
// Child class: Resistor, Inductor, Capacitor, Resistor_non_lin
class ConductanceDevice {
public:
    ConductanceDevice(Node* n1, Node* n2) : node1(n1), node2(n2) { }
    virtual std::complex<double> get_conductance(double omega) const = 0;

    Node* get_node1() const {
        return node1;
    }

    Node* get_node2() const {
        return node2;
    }

    virtual ~ConductanceDevice() { }

protected:
    Node* node1;
    Node* node2;

};

class Resistor : public ConductanceDevice {
public:

    Resistor(double r, Node* n1, Node* n2) : resistance(r), ConductanceDevice(n1, n2) { }

    std::complex<double> get_conductance(double omega) const {
        std::complex<double> conductance(1 / resistance);
        return conductance;
    }

private:
    double resistance;
};

class Capacitor : public ConductanceDevice {
public:
    Capacitor(double c, Node* n1, Node* n2) : capacitance(c), ConductanceDevice(n1, n2) { }

    std::complex<double>get_conductance(double omega) const {
        std::complex<double> conductance(0, omega * capacitance);
        return conductance;
    }

private:
    double capacitance;
};

class Inductor : public ConductanceDevice {
public:
    Inductor(double l, Node* n1, Node* n2) : inductance(l), ConductanceDevice(n1, n2) { }

    std::complex<double>get_conductance(double omega) const {
        std::complex<double> conductance(0, -1 / (omega * inductance));
        return conductance;
    }

private:
    double inductance;
};


// Inheritance chain: ConductanceDevices -> Resistor_non_lin_equiv -> Resistor_D & Resistor_Qbe & Resistor_Qce & Resistor Qbc & Resistor_M
class Resistor_non_lin_equiv : public ConductanceDevice {
public:

    Resistor_non_lin_equiv(Node* n1, Node* n2) : ConductanceDevice(n1, n2) { }
    virtual Node* get_n_a_b_g() const = 0;
    virtual Node* get_n_0_e_s() const = 0;
    virtual Node* get_n_c_c_d() const = 0;
    virtual void find_conductance(std::complex<double> v_a_b_g, std::complex<double> v_0_e_s, std::complex<double> v_c_c_d) = 0;
    std::complex<double> get_conductance(double omega) const {
        return conductance;
    }
protected:
    std::complex<double> conductance;
};

class Resistor_D : public Resistor_non_lin_equiv {
public:
    Resistor_D(Node* n1, Node* n2, Diode* diode) : Resistor_non_lin_equiv(n1, n2), D(diode) { }

    void find_conductance(std::complex<double> v_anode, std::complex<double> v_not_used, std::complex<double> v_cathode) {
        std::string mode = mode_diode(v_anode, v_cathode);
        if (mode == "on") {
            conductance = (Is / Vt) * exp((v_anode - v_cathode) / Vt);

        }
        else if (mode == "off") {
            conductance = 0;
        }
    }

    Node* get_n_a_b_g() const { // anode, base or gate
        return D->get_anode();
    }
    Node* get_n_0_e_s() const { // not-used, emitter or source
        return NULL;
    }
    Node* get_n_c_c_d() const { // cathode, collector or drain
        return D->get_cathode();
    }

private:
    Diode* D;
};

class Resistor_Qbe : public Resistor_non_lin_equiv {
public:
    Resistor_Qbe(Node* n1, Node* n2, BJT* bjt) : Resistor_non_lin_equiv(n1, n2), Q(bjt) { }

    void find_conductance(std::complex<double> vb, std::complex<double> ve, std::complex<double> vc) {
        std::string type = Q->get_model();
        if (type == "NPN") {
            std::string mode = mode_bjt(type, vc, ve, vb);
            if (mode == "active") {
                conductance = (Is / (Q->get_beta_forward() * Vt)) * exp((vb - ve) / Vt) * (1.0 + ((vc - ve) / Q->get_VA()));
            }
            else if (mode == "off") {
                conductance = 0;
            }
            else if (mode == "saturation") {
                conductance = Is / (Q->get_beta_forward() * Vt) * exp((vb - ve) / Vt);
            }
            else if (mode == "reverse") {
                conductance = 0;
            }
        }
        else if (type == "PNP") {
            std::string mode = mode_bjt(type, vc, ve, vb);
            if (mode == "active") {
                conductance = (Is / (Q->get_beta_forward() * Vt)) * exp((ve - vb) / Vt) * (1.0 + ((ve - vc) / Q->get_VA()));
            }
        }
    }

    Node* get_n_a_b_g() const { // anode, base or gate
        return Q->get_base();
    }
    Node* get_n_0_e_s() const { // not-used, emitter or source
        return Q->get_emitter();
    }
    Node* get_n_c_c_d() const { // cathode, collector or drain
        return Q->get_collector();
    }

private:
    BJT* Q;
};

class Resistor_Qce : public Resistor_non_lin_equiv {
public:
    Resistor_Qce(Node* n1, Node* n2, BJT* bjt) : Resistor_non_lin_equiv(n1, n2), Q(bjt) { }

    void find_conductance(std::complex<double> vb, std::complex<double> ve, std::complex<double> vc) {
        std::string type = Q->get_model();
        if (type == "NPN") {
            std::string mode = mode_bjt(type, vc, ve, vb);
            if (mode == "active") {
                conductance = Is * exp((vb - ve) / Vt) / Q->get_VA();
            }
            else if (mode == "off") {
                conductance = 0;
            }
            else if (mode == "saturation") {
                conductance = 0;
            }
            else if (mode == "reverse") {
                conductance = 0;
            }
        }
        else if (type == "PNP") {
            std::string mode = mode_bjt(type, vc, ve, vb);
            if (mode == "active") {
                conductance = Is * exp((ve - vb) / Vt) / Q->get_VA();
            }
        }

    }

    Node* get_n_a_b_g() const { // anode, base or gate
        return Q->get_base();
    }
    Node* get_n_0_e_s() const { // not-used, emitter or source
        return Q->get_emitter();
    }
    Node* get_n_c_c_d() const { // cathode, collector or drain
        return Q->get_collector();
    }

private:
    BJT* Q;
};

class Resistor_Qbc : public Resistor_non_lin_equiv {
public:
    Resistor_Qbc(Node* n1, Node* n2, BJT* bjt) : Resistor_non_lin_equiv(n1, n2), Q(bjt) { }

    void find_conductance(std::complex<double> vb, std::complex<double> ve, std::complex<double> vc) {
        std::string type = Q->get_model();
        if (type == "NPN") {
            std::string mode = mode_bjt(type, vc, ve, vb);
            double betar = Q->get_beta_reverse();
            double alphar = betar / (betar + 1);
            if (mode == "active") {
                conductance = 0;
            }
            else if (mode == "off") {
                conductance = 0;
            }
            else if (mode == "saturation") {
                conductance = (Is / (Vt * alphar)) * exp((vb - vc) / Vt);
            }
            else if (mode == "reverse") {
                conductance = (Is / (betar * Vt)) * exp((vb - vc) / Vt);
            }
        }
        else if (type == "PNP") {
            std::string mode = mode_bjt(type, vc, ve, vb);
            if (mode == "active") {
                conductance = 0;
            }
        }
    }

    Node* get_n_a_b_g() const { // anode, base or gate
        return Q->get_base();
    }
    Node* get_n_0_e_s() const { // not-used, emitter or source
        return Q->get_emitter();
    }
    Node* get_n_c_c_d() const { // cathode, collector or drain
        return Q->get_collector();
    }

private:
    BJT* Q;
};

class Resistor_M : public Resistor_non_lin_equiv {
public:
    Resistor_M(Node* n1, Node* n2, Mosfet* mosfet) : Resistor_non_lin_equiv(n1, n2), M(mosfet) { }

    void find_conductance(std::complex<double> vg, std::complex<double> vs, std::complex<double> vd) {
        std::string type = M->get_model();
        std::string mode = mode_mosfet(type, vd, vg, vs, M->get_Vt_mosfet());
        if (mode == "linear") {
            conductance = 2 * M->get_k() * ((vg - vs) - M->get_Vt_mosfet() - (vd - vs));

        }
        else if (mode == "saturation") {
            conductance = M->get_k() * pow((vg - vs - (M->get_Vt_mosfet())), 2) / M->get_VA();
        }
        else if (mode == "off") {
            conductance = 0;
        }
    }
    Node* get_n_a_b_g() const { // anode, base or gate
        return M->get_gate();
    }
    Node* get_n_0_e_s() const { // not-used, emitter or source
        return M->get_source();
    }
    Node* get_n_c_c_d() const { // cathode, collector or drain
        return M->get_drain();
    }
private:
    Mosfet* M;
};


class Vcontrolled_Isource {
public:
    //current flows from n_pos to n_neg;
    //n_a_b_g is anode/base/gate
    //n_0_e_s is 0/emitter/source
    //n_c_c_d is cathode/collector/drain
    Vcontrolled_Isource(Node* n_pos, Node* n_neg, Node* n_con_pos, Node* n_con_neg) : pos(n_pos), neg(n_neg), con_pos(n_con_pos), con_neg(n_con_neg) { }
    virtual void find_gm(std::complex<double> v_a_b_g, std::complex<double> v_0_e_s, std::complex<double> v_c_c_d) = 0;
    virtual Node* get_n_a_b_g() const = 0;
    virtual Node* get_n_0_e_s() const = 0;
    virtual Node* get_n_c_c_d() const = 0;
    Node* get_neg() const {
        return neg;
    }
    Node* get_pos() const {
        return pos;
    }
    Node* get_con_neg() const {
        return con_neg;
    }
    Node* get_con_pos() const {
        return con_pos;
    }
    std::complex<double> get_gm() const {
        return gm;
    }

protected:
    Node* pos;
    Node* neg;
    Node* con_pos;
    Node* con_neg;
    std::complex<double> gm;
};

//Voltage controlled current source connected to CE (C is pos, E is neg) controlled by voltage across BE(B is con_pos E is con_neg)
class G_bjt_be : public Vcontrolled_Isource {
public:
    G_bjt_be(Node* n_pos, Node* n_neg, Node* n_con_pos, Node* n_con_neg, BJT* bjt) : Vcontrolled_Isource(n_pos, n_neg, n_con_pos, n_con_neg), Q(bjt) {}

    void find_gm(std::complex<double> v_b, std::complex<double> v_e, std::complex<double> v_c) {
        std::string type = Q->get_model();

        if (type == "NPN") {
            std::string mode = mode_bjt(type, v_c, v_e, v_b);
            if (mode == "active") {
                // active mode
                gm = (Is / Vt) * exp((v_b - v_e) / Vt) * (1.0 + (v_c - v_e) / Q->get_VA());
                
            }
            else if (mode == "saturation") {
                // saturation mode
                gm = (Is / Vt) * exp((v_b - v_e) / Vt);
            }
            else if (mode == "off") {
                // off mode
                gm = 0;
            }
            else if (mode == "reverse") {
                // reverse mode
                gm = 0;
            }
        }
        else if (type == "PNP") {
            std::string mode = mode_bjt(type, v_c, v_e, v_b);
            if (mode == "active") {
                // active mode
                gm = (Is / Vt) * exp((v_e - v_b) / Vt) * (1.0 + (v_e - v_c) / Q->get_VA());
                
            }
        }
    }
    Node* get_n_a_b_g() const { // anode, base or gate
        return Q->get_base();
    }

    Node* get_n_0_e_s() const { // not-used, emitter or source
        return Q->get_emitter();
    }
    Node* get_n_c_c_d() const { // cathode, collector or drain
        return Q->get_collector();
    }

private:
    BJT* Q;
};

//Voltage controlled current source connected to EC (E is pos,  Cis neg) controlled by voltage across BC (B is con_pos C is con_neg)
class G_bjt_bc : public Vcontrolled_Isource {
public:
    G_bjt_bc(Node* n_pos, Node* n_neg, Node* n_con_pos, Node* n_con_neg, BJT* bjt) : Vcontrolled_Isource(n_pos, n_neg, n_con_pos, n_con_neg), Q(bjt) {}
    void find_gm(std::complex<double> v_b, std::complex<double> v_e, std::complex<double> v_c) {
        std::string type = Q->get_model();
        if (type == "NPN") {
            std::string mode = mode_bjt(type, v_c, v_e, v_b);
            if (mode == "active") {
                // active mode
                gm = 0;
            }
            else if (mode == "saturation") {
                // saturation mode
                gm = 0;
            }
            else if (mode == "off") {
                // off mode
                gm = 0;
            }
            else if (mode == "reverse") {
                // reverse mode
                gm = (Is / Vt) * exp((v_b - v_c) / Vt);
            }
        }
        else if (type == "PNP") {
            std::string mode = mode_bjt(type, v_c, v_e, v_b);
            if (mode == "active") {
                gm = 0;
            }
        }

    }
    Node* get_n_a_b_g() const { // anode, base or gate
        return Q->get_base();
    }

    Node* get_n_0_e_s() const { // not-used, emitter or source
        return Q->get_emitter();
    }
    Node* get_n_c_c_d() const { // cathode, collector or drain
        return Q->get_collector();
    }

private:
    BJT* Q;
};

//Voltage controlled current source connected to DS (D is pos, S is neg) controlled by voltage across GS (G is con_pos, S is con_neg)
class G_mosfet : public Vcontrolled_Isource {
public:
    G_mosfet(Node* n_pos, Node* n_neg, Node* n_con_pos, Node* n_con_neg, Mosfet* mosfet) : Vcontrolled_Isource(n_pos, n_neg, n_con_pos, n_con_neg), M(mosfet) {}
    void find_gm(std::complex<double> v_g, std::complex<double> v_s, std::complex<double> v_d) {
        std::complex<double> gm;
        std::string type = M->get_model();
        std::string mode = mode_mosfet(type, v_d, v_g, v_s, M->get_Vt_mosfet());
        if (mode == "linear") {
            // linear/triode mode
            gm = 2 * M->get_k() * (v_d - v_s);
        }
        else if (mode == "saturation") {
            // saturation mode
            gm = 2 * M->get_k() * (v_g - v_s - M->get_Vt_mosfet());
        }
        else if (mode == "off") {
            // off mode
            gm = 0;
        }
    }
    Node* get_n_a_b_g() const { // anode, base or gate
        return M->get_gate();
    }
    Node* get_n_0_e_s() const { // not-used, emitter or source
        return M->get_source();
    }
    Node* get_n_c_c_d() const { // cathode, collector or drain
        return M->get_drain();
    }

private:
    Mosfet* M;
};


// Parent class: Source
// Child class: Voltage_source, Current_source
class Source {
public:

    Source(Node* n1, Node* n2) : pos(n1), neg(n2) { }

    ~Source();

    Node* get_n_pos() const {
        return pos;
    }

    Node* get_n_neg() const {
        return neg;
    }

protected:
    Node* pos;
    Node* neg;
};


class Voltage_source : public Source {
public:
    Voltage_source(Node* n1, Node* n2, std::complex<double> voltage) : Source(n1, n2), V(voltage) {
    }

    std::complex<double> get_V() const {
        return V;
    }

private:
    std::complex<double> V;
};


class Current_source : public Source {
public:
    Current_source(Node* n1, Node* n2, std::complex<double> current) : Source(n1, n2), I(current) {
    }
    void change_I(int new_I) { I = new_I; }
    std::complex<double> get_I() const {
        return I;
    }

private:
    std::complex<double> I;
};


class iSource_nonlin {
public:

    iSource_nonlin(Node* n_pos, Node* n_neg) : pos(n_pos), neg(n_neg) { }

    Node* const get_pos() const {
        return pos;
    }

    Node* const get_neg() const {
        return neg;
    }

    std::complex<double> get_I() const {
        return I;
    }
    virtual void find_I(std::complex<double> v_a_b_g, std::complex<double> v_0_e_s, std::complex<double> v_c_c_d) = 0;
    virtual Node* get_n_a_b_g() const = 0;
    virtual Node* get_n_0_e_s() const = 0;
    virtual Node* get_n_c_c_d() const = 0;

    ~iSource_nonlin();

protected:
    Node* pos;
    Node* neg;
    std::complex<double> I;
};

//flows from anode to cathode
class I_source_D : public iSource_nonlin {
public:
    I_source_D(Node* pos, Node* neg) : iSource_nonlin(pos, neg) { }
    void find_I(std::complex<double> v_a, std::complex<double> v_0, std::complex<double> v_c) {
        std::string mode = mode_diode(v_a, v_c);
        if (mode == "on") {
            I = Is * (exp((v_a - v_c) / Vt) - 1.0) - (Is / Vt) * exp((v_a - v_c) / Vt) * (v_a - v_c);
        }
        else if (mode == "off") {
            I = 0;
        }
    }
    Node* get_n_a_b_g() const {
        return pos;
    }
    Node* get_n_0_e_s() const {
        return NULL;
    }
    Node* get_n_c_c_d() const {
        return neg;
    }
};

//flows from drain to source, drain is pos, source is neg
class I_source_mosfet :public iSource_nonlin {
public:
    I_source_mosfet(Node* n_pos, Node* n_neg, Mosfet* mosfet) : iSource_nonlin(n_pos, n_neg), M(mosfet) { }
    void find_I(std::complex<double> v_g, std::complex<double> v_s, std::complex<double> v_d) {
        std::string type = M->get_model();
        double VT = M->get_Vt_mosfet();
        std::string mode = mode_mosfet(type,v_d, v_g, v_s, VT);
        double K = M->get_k();
        double VA = M->get_VA();
        if (mode == "linear") {
            I = K * (2.0 * (v_g - v_s - VT) * (v_d - v_s) - pow((v_d - v_s), 2)) - 2 * K * (v_d - v_s) * (v_g - v_s) - 2.0 * K * (v_g - v_s - VT - v_d + v_s) * (v_d - v_s);

        }
        else if (mode == "saturation") {
            I = K * pow(v_g - v_s - VT, 2) * (1.0 + ((v_d - v_s) / VA)) - 2 * K * (v_g - v_s - VT) * (v_g - v_s) - (K / VA) * std::pow((v_g - v_s - VT), 2) * (v_d - v_s);
            
        }
        else if (mode == "off") {
            I = 0;
        }

    }
    Node* get_n_a_b_g() const { // anode, base or gate
        return M->get_gate();
    }
    Node* get_n_0_e_s() const { // not-used, emitter or source
        return M->get_source();
    }
    Node* get_n_c_c_d() const { // cathode, collector or drain
        return M->get_drain();
    }

private:
    Mosfet* M;
};

//from C to E (C is pos E is neg)
class I_source_Qce :public iSource_nonlin {
public:
    I_source_Qce(Node* n_pos, Node* n_neg, BJT* bjt) :iSource_nonlin(n_pos, n_neg), Q(bjt) { }
    void find_I(std::complex<double> v_b, std::complex<double> v_e, std::complex<double> v_c) {
        std::string type = Q->get_model();
        std::complex<double> Vce = v_c - v_e;
        std::complex<double> Vec = v_e - v_c;
        std::complex<double> Vbe = v_b - v_e;
        std::complex<double> Veb = v_e - v_b;
        std::complex<double> Vbc = v_b - v_c;
        std::complex<double> Vcb = v_c - v_b;
        double VA = Q->get_VA();
        double betar = Q->get_beta_reverse();
        double alphar = betar / (betar + 1);
        if (type == "NPN") {
            std::string mode = mode_bjt(type, v_c, v_e, v_b);

            if (mode == "active") {
                I = Is * exp(Vbe / Vt) * (1.0 + Vce / VA) - (Is / Vt) * exp(Vbe / Vt) * (1.0 + Vce / VA) * Vbe - Is * (exp(Vbe / Vt) / VA) * Vce;
            }
            else if (mode == "saturation") {

                I = Is * exp(Vbe / Vt) - (Is / alphar) * exp(Vbc / Vt) - (Is / Vt) * exp(Vbe / Vt) * Vbe - (Is / (Vt * alphar)) * exp(Vbc / Vt) * Vcb;
            }
            else if (mode == "off") {
                I = 0;
            }
            else if (mode == "reverse") {
                I = (Is / Vt) * exp(Vbc / Vt) * Vbc - Is * exp(Vbc / Vt);

            }
        }
        else if (type == "PNP") {
            std::string mode = mode_bjt(type, v_c, v_e, v_b);
            if (mode == "active") {
                I = Is * exp(Veb / Vt) * (1.0 + Vec / VA) - (Is / Vt) * exp(Veb / Vt) * (1.0 + Vec / VA) * Veb - Is * (exp(Veb / Vt) / VA) * Vec;
            }
        }


    }
    Node* get_n_a_b_g() const { // anode, base or gate
        return Q->get_base();
    }
    Node* get_n_0_e_s() const { // not-used, emitter or source
        return Q->get_emitter();
    }
    Node* get_n_c_c_d() const { // cathode, collector or drain
        return Q->get_collector();
    }

private:
    BJT* Q;
};

//flow from B to E (i.e. B is pos, E is neg)
class I_source_Qbe :public iSource_nonlin {
public:
    I_source_Qbe(Node* n_pos, Node* n_neg, BJT* bjt) :iSource_nonlin(n_pos, n_neg), Q(bjt) { }
    void find_I(std::complex<double> v_b, std::complex<double> v_e, std::complex<double> v_c) {
        std::complex<double> Vce = v_c - v_e;
        std::complex<double> Vec = v_e - v_c;
        std::complex<double> Vbe = v_b - v_e;
        std::complex<double> Veb = v_e - v_b;
        std::complex<double> Vbc = v_b - v_c;
        double VA = Q->get_VA();
        double betaf = Q->get_beta_forward();
        double betar = Q->get_beta_reverse();
        double alphar = betar / (betar + 1);
        std::string type = Q->get_model();
        if (type == "NPN") {
            std::string mode = mode_bjt(type, v_c, v_e, v_b);

            if (mode == "active") {
                I = Is * exp(Vbe / Vt) * (1.0 + Vce / VA) / betaf - (1.0 + Vce / VA) * (Is / Vt) * exp(Vbe / Vt) * (Vbe / betaf);
      
            }
            else if (mode == "saturation") {

                I = (Is / betaf) * exp(Vbe / Vt) + (Is / alphar) * exp(Vbc / Vt) - (Is / (Vt * betaf)) * exp(Vbe / Vt) * Vbe - (Is / (alphar * Vt)) * exp(Vbc / Vt) * Vbc;
              
            }
            else if (mode == "off") {
                I = 0;
            }
            else if (mode == "reverse") {
                I = 0;
            }
        }
        else if (type == "PNP") {
            std::string mode = mode_bjt(type, v_c, v_e, v_b);
            if (mode == "active") {
                I = Is * exp(Veb / Vt) * (1.0 + Vec / VA) / betaf - (1.0 + Vec / VA) * (Is / Vt) * exp(Veb / Vt) * (Veb / betaf);
               
            }
        }

    }
    Node* get_n_a_b_g() const { // anode, base or gate
        return Q->get_base();
    }
    Node* get_n_0_e_s() const { // not-used, emitter or source
        return Q->get_emitter();
    }
    Node* get_n_c_c_d() const { // cathode, collector or drain
        return Q->get_collector();
    }

private:
    BJT* Q;
};

//flows from B to C (B is pos, C is neg)
class I_source_Qbc :public iSource_nonlin {
public:
    I_source_Qbc(Node* n_pos, Node* n_neg, BJT* bjt) :iSource_nonlin(n_pos, n_neg), Q(bjt) { }
    void find_I(std::complex<double> v_b, std::complex<double> v_e, std::complex<double> v_c) {
        std::string type = Q->get_model();
        if (type == "NPN") {
            std::string mode = mode_bjt(type, v_c, v_e, v_b);
            std::complex<double> Vbc = v_b - v_c;
            double betar = Q->get_beta_reverse();
            if (mode == "active") {
                I = 0;
            }
            else if (mode == "saturation") {
                I = 0;
            }
            else if (mode == "off") {
                I = 0;
            }
            else if (mode == "reverse") {
                I = (Is / betar) * (exp(Vbc / Vt) - 1.0) - (Is / (betar * Vt)) * exp(Vbc / Vt) * Vbc;
            }
        }
        else if (type == "PNP") {
            std::string mode = mode_bjt(type, v_c, v_e, v_b);
            if (mode == "active") {
                I = 0;
            }
        }
    }
    Node* get_n_a_b_g() const { // anode, base or gate
        return Q->get_base();
    }
    Node* get_n_0_e_s() const { // not-used, emitter or source
        return Q->get_emitter();
    }
    Node* get_n_c_c_d() const { // cathode, collector or drain
        return Q->get_collector();
    }

private:
    BJT* Q;
};

#endif










