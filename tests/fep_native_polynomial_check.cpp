// Optional native regression driver; compile with the same core as the extension.
#include "core.hpp"
#include <cmath>
#include <iostream>
using namespace xpongecpp;

Molecule make_molecule() {
    Molecule m("native-regression");
    Residue residue; residue.name="NATIVE";residue.atom_count=4;m.residues.push_back(residue);
    for (int i=0;i<4;++i) {Atom a;a.name="C"+std::to_string(i);a.type="C";a.element="C";a.mass=12.;a.x=i;a.y=i%2;a.z=.1*i;m.atoms.push_back(a);}
    m.topology_override=Topology{};m.topology_override->exclusions.resize(4);
    return m;
}

double bond_energy(const Topology& topology,double r) {
    double out=0.;for (const auto& f:topology.bonds)out+=f.k*(r-f.length)*(r-f.length);return out;
}
double angle_energy(const Topology& topology,double theta) {
    double out=0.;for (const auto& f:topology.angles)out+=f.k*(theta-f.theta)*(theta-f.theta);return out;
}
double dihedral_energy(const Topology& topology,double phi) {
    double out=0.;for (const auto& f:topology.dihedrals)out+=f.k*(1+std::cos(f.periodicity*phi-f.phase));return out;
}

int main() {
    int failed=0,total=0;
    const auto check=[&](const char* name,bool passed) {++total;failed+=!passed;std::cout<<name<<": "<<(passed?"PASS":"FAIL")<<'\n';};
    auto a=make_molecule(),b=a;
    a.topology_override->bonds={{0,1,100.,1.3},{0,1,200.,1.7},{0,1,50.,1.3}};
    b.topology_override->bonds={{1,0,300.,1.5},{0,1,80.,1.3}};
    a.topology_override->angles={{0,1,2,20.,1.2},{0,1,2,30.,1.8}};
    b.topology_override->angles={{2,1,0,40.,1.4},{0,1,2,50.,1.8}};
    a.topology_override->dihedrals={{0,1,2,3,2,2.,.1234561},{0,1,2,3,2,3.,.1234562},{0,1,2,3,2,7.,.1234561}};
    b.topology_override->dihedrals={{3,2,1,0,2,8.,.1234562},{0,1,2,3,3,4.,.7}};
    bool bond_ok=true,angle_ok=true,dihedral_ok=true,derivative_ok=true;
    for (double l : {0.,.2,.5,.8,1.}) {
        const auto merged=build_topology(merge_force_field(a,b,l,{}));
        for (double x : {-.7,.2,1.4,2.1,3.}) {
            bond_ok &= std::abs(bond_energy(merged,x)-(1-l)*bond_energy(*a.topology_override,x)-l*bond_energy(*b.topology_override,x))<1e-9;
            angle_ok &= std::abs(angle_energy(merged,x)-(1-l)*angle_energy(*a.topology_override,x)-l*angle_energy(*b.topology_override,x))<1e-9;
            dihedral_ok &= std::abs(dihedral_energy(merged,x)-(1-l)*dihedral_energy(*a.topology_override,x)-l*dihedral_energy(*b.topology_override,x))<1e-9;
            const double h=1e-5;
            const auto d=[&](const Topology& t){return (dihedral_energy(t,x+h)-dihedral_energy(t,x-h))/(2*h);};
            derivative_ok &= std::abs(d(merged)-(1-l)*d(*a.topology_override)-l*d(*b.topology_override))<1e-8;
        }
    }
    check("multiple_bond_equilibria_and_duplicates",bond_ok);
    check("multiple_angle_equilibria_and_duplicates",angle_ok);
    check("close_phases_and_duplicate_fourier_terms",dihedral_ok);
    check("finite_difference_of_torsion_function",derivative_ok);
    const auto middle=build_topology(merge_force_field(a,b,.5,{}));
    check("close_phases_remain_distinct",middle.dihedrals.size()==3);
    auto c=make_molecule(),d=c;c.topology_override->dihedrals={{0,1,2,3,2,2.,.2}};
    d.topology_override->dihedrals={{3,2,1,0,2,6.,.2}};
    const auto reverse=build_topology(merge_force_field(c,d,.25,{}));
    check("full_reversal_is_one_equivalent_term",reverse.dihedrals.size()==1&&std::abs(reverse.dihedrals[0].k-3.)<1e-12);
    std::cout<<"total="<<total<<" failed="<<failed<<'\n';return failed?1:0;
}
