#ifndef short_circuit_h
#define short_circuit_h

void short_v_source_dc(std::vector<Voltage_source*> vsources, std::vector<Node*>& nodes);
void short_inductor(std::vector<Inductor*> inductors, std::vector<Node*>& nodes);
void short_v_source_ac(std::vector<Voltage_source*> vsources, std::vector<Node*>& nodes);

#endif


