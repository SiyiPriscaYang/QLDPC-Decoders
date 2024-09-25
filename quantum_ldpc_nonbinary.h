//
// Created by Siyi Yang on 2/15/24.
//

#ifndef QLDPC_DEGENERATE_QUANTUM_LDPC_NONBINARY_H
#define QLDPC_DEGENERATE_QUANTUM_LDPC_NONBINARY_H

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <string>
#include "tools.h"
#include "node_basic.h"


class qVN: public Node{
protected:
    vector<int> stabilizer;
    vector<vector<double> >  qmessage;
public:
    qVN(){vector<int> neigh(0);vector<double> mes(0);vector<int> stab(0);vector<vector<double> >  qmes(0);neighbors=neigh;message=mes;stabilizer=stab;qmessage=qmes;}
    qVN(vector<int> neigh, vector<double> mes,vector<int> stab,vector<vector<double> >  qmes){neighbors=neigh;message=mes;stabilizer=stab;qmessage=qmes;}
    void set_qmessage(int ind,double x,double y,double z){qmessage[ind][0]=x;qmessage[ind][1]=y;qmessage[ind][2]=z;}
    vector<double> update_v2c(vector<double> lh);// update v2c messages
    void print2(){
        for(int i=0;i<qmessage.size();i++) {
            cout<<i<<':';
            for (int j = 0; j < 3; j++)
                cout << qmessage[i][j] << ',';
            cout<<endl;
        }
    }
};

class qCode{
private:
    vector<qVN> vn;// store VNs
    vector<CN> cn;// store CNs
    vector<vector<int> >  vc;//store VN indices adjacent to CNs
    vector<vector<int> >  sc;//store VN stabilizer values adjacent to CNs
    vector<double> lv2c;// store v2c messages of all edges
    vector<double> lc2v;// store c2v messages of all edges
    vector<vector<int> > bx;
    vector<vector<int> > bz;
public:
    qCode(vector<vector<int> >  H);
    void setxz(vector<vector<int> > H);
    vector<int> get_syndrome(vector<int> x);
    bool is_satisfy(vector<int> x,vector<int> syndrome);
    void update_c2v(vector<int> syndrome);
    vector<vector<double> >  update_v2c(vector<vector<double> >  lch);
    vector<int> decode(vector<int> x, vector<vector<double> >  lch, int num_iter);
    bool degenerate(vector<int> x);
    void print(){
        for(int i=0;i<cn.size();i++) {
            cout<<"cn-"<<i<<':';
            cn[i].print();
        }
        for(int i=0;i<vn.size();i++) {
            cout<<"vn-"<<i<<':';
            vn[i].print();
        }
    }
};


// checked
vector<double> qVN::update_v2c(vector<double> lh) {
    vector<double> sum(lh);
    int i,j;
    double a;
    for(i=0;i<message.size();i++) {
        a=message[i];
        if(stabilizer[i]==1)
            set_qmessage(i,0,-a,-a);// a=log(p0/p1)=log(pi+px/py+pz)=log(pi/py,pz)=log(px/py,pz)
        if(stabilizer[i]==2)
            set_qmessage(i,-a,0,-a);
        if(stabilizer[i]==3)
            set_qmessage(i,-a,-a,0);
        for (j = 0; j < 3; j++)
            sum[j] += qmessage[i][j];
        //cout<<"stabilizer="<<stabilizer[i]<<','<<a<<':';
        //print_double(sum);
    }

    for(i=0;i<message.size();i++) {
        for (j = 0; j < 3; j++)
            qmessage[i][j] = sum[j] - qmessage[i][j];
        if(stabilizer[i]==1)
            a=q2b(0,qmessage[i][0])-q2b(qmessage[i][1],qmessage[i][2]);
        if(stabilizer[i]==2)
            a=q2b(0,qmessage[i][1])-q2b(qmessage[i][2],qmessage[i][0]);
        if(stabilizer[i]==3)
            a=q2b(0,qmessage[i][2])-q2b(qmessage[i][0],qmessage[i][1]);
        set_message(i,a);
    }
    return sum;
}


void qCode::setxz(vector<vector<int> > H) {
    int m=H.size();
    int n=H[0].size();
    int t;
    vector<vector<int> > Hx(m, vector<int>(n, 0));
    vector<vector<int> > Hz(m, vector<int>(n, 0));
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            t=H[i][j];
            if(t){
                Hx[i][j]=((t==1)||(t==2))? 1:0;
                Hz[i][j]=((t==2)||(t==3))? 1:0;
            }
        }
    }

    bx= generator(Hx);
    bz= generator(Hz);
}

bool qCode::degenerate(vector<int> x) {
    int len=x.size();
    vector<int> vx(len);
    vector<int> vz(len);
    int sys=bx.size();
    vector<int> sysx(0);
    vector<int> sysz(0);
    int t;

    for(int i=0;i<len;i++){
        t=x[i];
        vx[i]=((t==1)||(t==2))? 1:0;
        vz[i]=((t==2)||(t==3))? 1:0;
    }

    vector<int> idx= is_in_row_space(bx,vx);
    if(idx.size())
        if(idx[0]==-1)
            return false;
    vector<int> idz= is_in_row_space(bz,vz);
    if(idz.size())
        if(idz[0]==-1)
            return false;
/*
    vector<int> idx1, idz1;
    set_difference(idx.begin(), idx.end(), idz.begin(), idz.end(), back_inserter(idx1));
    set_difference(idz.begin(), idz.end(), idx.begin(), idx.end(), back_inserter(idz1));

    for(int i: idx1){
        for(int j=0;j<len;j++){
            if(bx[i][j]&&(!bz[i][j]))
                return false;
        }
    }
    for(int i: idz1){
        for(int j=0;j<len;j++){
            if(bz[i][j]&&(!bx[i][j]))
                return false;
        }
    }
*/

    return true;
}


qCode::qCode(vector<vector<int> >  H){
    int num_cn=H.size();
    int num_vn=H[0].size();
    int num_e=0;
    int e=0;
    vector<vector<int> >  Ht=H;
    vector<vector<int> >  vct(num_cn);
    vector<vector<int> >  sct(num_cn);
    vector<qVN> vt(num_vn);
    vector<CN> ct(num_cn);
    for(int c=0;c<num_cn;c++){
        vector<int> neigh(0);// neigh_cn in e
        vector<int> neigh_c(0);// neigh_cn in v
        vector<int> neigh_sc(0);// neigh_cn in stabilizer
        for(int v=0;v<num_vn;v++)
            if(H[c][v]) {
                neigh.push_back(num_e);// e
                neigh_c.push_back(v);
                neigh_sc.push_back(H[c][v]);
                Ht[c][v]=num_e;// Ht records edge number
                num_e++;
                //cout<<H[c][v]<<',';
            }
        //cout<<endl;
        vct[c]=neigh_c;// in v
        sct[c]=neigh_sc;// in stabilizer value
        vector<double> mes(neigh.size());
        CN node(neigh,mes);// check
        ct[c]=node;
    }

    for(int v=0;v<num_vn;v++){
        vector<int> neigh(0);
        vector<int> stab(0);
        for(int c=0;c<num_cn;c++)
            if(H[c][v]) {
                e=Ht[c][v];
                neigh.push_back(e);
                stab.push_back(H[c][v]);
            }
        vector<double> mes(neigh.size());
        vector<vector<double> >  qmes(neigh.size());
        vector<double> qmest(3);
        for(int i=0;i<neigh.size();i++)
            qmes[i]=qmest;
        qVN node(neigh,mes,stab,qmes);
        vt[v]=node;
    }

    vector<double> c2v(num_e);
    vector<double> v2c(num_e);

    vn=vt;
    cn=ct;
    lv2c=v2c;
    lc2v=c2v;
    vc=vct;
    sc=sct;
    setxz(H);
}


// checked
vector<int> qCode::get_syndrome(vector<int> x) {
    vector<int> syndrome(vc.size());
    vector<int> lv,ls;
    int parity,ind,t,s;
    for(int c=0;c<vc.size();c++){
        parity=0;
        lv=vc[c];
        ls=sc[c];
        for(int v=0;v<lv.size();v++) {
            ind=lv[v];
            s=ls[v];// edge pauli
            t=x[ind];// vn pauli
            if(t&&(t!=s))
                parity = parity ^ 1;
            //if(t)
            //    cout<<v<<':'<<s<<','<<x[ind]<<';'<<"parity="<<parity<<',';
        }
        syndrome[c]=parity;
        //if(syndrome[c])
        //    cout<<"syndrome["<<c<<"]="<<parity<<';';
    }
    //cout<<endl;
    return syndrome;
}

// checked
bool qCode::is_satisfy(vector<int> x,vector<int> syndrome) {
    vector<int> lv,ls;
    int parity,ind,t,s;
    for(int c=0;c<vc.size();c++){
        parity=0;
        lv=vc[c];
        ls=sc[c];
        for(int v=0;v<lv.size();v++) {
            ind=lv[v];
            s=ls[v];
            t=x[ind];
            if(t&&(t!=s))
                parity = parity ^ 1;
        }
        if(parity!=syndrome[c])
            return false;
    }

    return true;
}



// checked
vector<vector<double> >  qCode::update_v2c(vector<vector<double> >  lch) {
    vector<double> tv2c;
    vector<vector<double> >  lh(vn.size());
    vector<int> lc;
    for(int n=0;n<vn.size();n++){
        lc=vn[n].get_neighbors();// vn neighbor in e
        for(int c=0;c<lc.size();c++)
            vn[n].set_message(c,lc2v[lc[c]]);
        lh[n]=vn[n].update_v2c(lch[n]);// already changed mes
        tv2c=vn[n].get_message();// should change to get qmes
        for(int c=0;c<lc.size();c++)
            lv2c[lc[c]]=tv2c[c];
    }
    return lh;
}



// checked
void qCode::update_c2v(vector<int> syndrome) {
    vector<double> tc2v;
    vector<int> lv;
    for(int n=0;n<cn.size();n++){
        lv=cn[n].get_neighbors();
        for(int v=0;v<lv.size();v++)
            cn[n].set_message(v,lv2c[lv[v]]);
        cn[n].update_c2v(syndrome[n]);
        tc2v=cn[n].get_message();
        for(int v=0;v<lv.size();v++)
            lc2v[lv[v]]=tc2v[v];
    }
}


// checked, err{{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}}
vector<int> qCode::decode(vector<int> x, vector<vector<double> >  lch, int num_iter) {
    vector<int> syndrome=get_syndrome(x);
    vector<int> diff(x.size());
    vector<vector<double> >  llr;
    vector<vector<int> > err{{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}};
    bool suc_decode;
    int decode_type;
    int l,v,t;
    bool satisfy;
    bool identical;
    vector<int> xt(x.size());

    for(l=0;l<num_iter;l++){
        llr=update_v2c(lch);
        update_c2v(syndrome);
        for(v=0;v<vn.size();v++)
            xt[v]=estimate(llr[v]);
        satisfy = is_satisfy(xt,syndrome);
        if(satisfy)
            break;
    }

    identical=true;
    for(v=0;v<vn.size();v++){
        if(xt[v]!=x[v]) {
            diff[v] = err[xt[v]][x[v]];
            if(identical)
                identical=false;
        }
        else
            diff[v]=0;
    }

    if(satisfy) {
        if (identical) {
            suc_decode=true;
            //cout << "Decoder succeed\n";
            diff.push_back(0);
        }
        else{
            suc_decode = degenerate(diff);
            if (suc_decode) {
                //cout << "Decoder succeed\n";
                diff.push_back(0);
            }
            else {
                //cout << "Decoder failed: convergent error\n";
                diff.push_back(1);// convergent error
            }
        }
    }
    else{
        suc_decode=false;
        //cout<<"Decoder failed: non-convergent error\n";
        diff.push_back(2);// non-convergent error
    }

    return diff;
}


#endif //QLDPC_DEGENERATE_QUANTUM_LDPC_NONBINARY_H
