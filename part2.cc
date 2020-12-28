#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <assert.h>
#include <math.h>
using namespace std;

void printUsage();
void genLayer(int N, int M, int T, int P, vector<int>& constvector, string modName, ofstream &os);
void genAllLayers(int N, int M1, int M2, int M3, int T, int A, vector<int>& constVector, string modName, ofstream &os);
void readConstants(ifstream &constStream, vector<int>& constvector);
void genROM(vector<int>& constVector, int bits, string modName, ofstream &os);




int main(int argc, char* argv[]) {

    // If the user runs the program without enough parameters, print a helpful message
    // and quit.
    if (argc < 7) {
        printUsage();
        return 1;
    }

    int mode = atoi(argv[1]);

    ifstream const_file;
    ofstream os;
    vector<int> constVector;

    //----------------------------------------------------------------------
    // Look here for Part 1 and 2
    if ((mode == 1) && (argc == 7)) {
        // Mode 1: Generate one layer with given dimensions and one testbench

        // --------------- read parameters, etc. ---------------
        int N = atoi(argv[2]);
        int M = atoi(argv[3]);
        int T = atoi(argv[4]);
        int P = atoi(argv[5]);
        const_file.open(argv[6]);
        if (const_file.is_open() != true) {
            cout << "ERROR reading constant file " << argv[6] << endl;
            return 1;
        }

        // Read the constants out of the provided file and place them in the constVector vector
        readConstants(const_file, constVector);

        string out_file = "conv_" + to_string(N) + "_" + to_string(M) + "_" + to_string(T) + "_" + to_string(P) + ".sv";

        os.open(out_file);
        if (os.is_open() != true) {
            cout << "ERROR opening " << out_file << " for write." << endl;
            return 1;
        }
        // -------------------------------------------------------------

        // call the genLayer function you will write to generate this layer
        string modName = "conv_" + to_string(N) + "_" + to_string(M) + "_" + to_string(T) + "_" + to_string(P);
        genLayer(N, M, T, P, constVector, modName, os);

    }
    //--------------------------------------------------------------------


    // ----------------------------------------------------------------
    // Look here for Part 3
    else if ((mode == 2) && (argc == 9)) {
        // Mode 2: Generate three layer with given dimensions and interconnect them

        // --------------- read parameters, etc. ---------------
        int N  = atoi(argv[2]);
        int M1 = atoi(argv[3]);
        int M2 = atoi(argv[4]);
        int M3 = atoi(argv[5]);
        int T  = atoi(argv[6]);
        int A  = atoi(argv[7]);
        const_file.open(argv[8]);
        if (const_file.is_open() != true) {
            cout << "ERROR reading constant file " << argv[8] << endl;
            return 1;
        }
        readConstants(const_file, constVector);

        string out_file = "multi_" + to_string(N) + "_" + to_string(M1) + "_" + to_string(M2) + "_" + to_string(M3) + "_" + to_string(T) + "_" + to_string(A) + ".sv";


        os.open(out_file);
        if (os.is_open() != true) {
            cout << "ERROR opening " << out_file << " for write." << endl;
            return 1;
        }
        // -------------------------------------------------------------

        string mod_name = "multi_" + to_string(N) + "_" + to_string(M1) + "_" + to_string(M2) + "_" + to_string(M3) + "_" + to_string(T) + "_" + to_string(A);

        // call the genAllLayers function
        genAllLayers(N, M1, M2, M3, T, A, constVector, mod_name, os);

    }
    //-------------------------------------------------------

    else {
        printUsage();
        return 1;
    }

    // close the output stream
    os.close();

}

// Read values from the constant file into the vector
void readConstants(ifstream &constStream, vector<int>& constvector) {
    string constLineString;
    while(getline(constStream, constLineString)) {
        int val = atoi(constLineString.c_str());
        constvector.push_back(val);
    }
}

// Generate a ROM based on values constVector.
// Values should each be "bits" number of bits.
void genROM(vector<int>& constVector, int bits, string modName, ofstream &os) {

        int numWords = constVector.size();
        int addrBits = ceil(log2(numWords));

        os << "module " << modName << "(clk, addr, z);" << endl;
        os << "   input                     clk;" << endl;
        os << "   input [" << addrBits-1 << ":0]               addr;" << endl;
        os << "   output logic signed [" << bits-1 << ":0] z;" << endl;
        os << "   always_ff @(posedge clk) begin" << endl;
        os << "      case(addr)" << endl;
        int i=0;
        for (vector<int>::iterator it = constVector.begin(); it < constVector.end(); it++, i++) {
            if (*it < 0)
                os << "        " << i << ": z <= -" << bits << "'d" << abs(*it) << ";" << endl;
            else
                os << "        " << i << ": z <= "  << bits << "'d" << *it      << ";" << endl;
        }
        os << "      endcase" << endl << "   end" << endl << "endmodule" << endl << endl;
}

// Parts 1 and 2
// Here is where you add your code to produce a neural network layer.
void genLayer(int N, int M, int T, int P, vector<int>& constVector, string modName, ofstream &os) {
    int addrbit_x = ceil(log2(N));
    int addrbit_f = ceil(log2(M));
    string romModName = modName + "_f_rom";
    int outbit = ceil(log2(N - M + 1));


    os<<"module " << modName << "(clk, reset, s_data_in_x, s_valid_x, s_ready_x, m_data_out_y, m_valid_y, m_ready_y);" << endl;
    os<<"   parameter T = 16;"<<endl;
    os<<"   input                      clk, reset, s_valid_x, m_ready_y;"<<endl;
    os<<"   input        signed ["<<T-1<<":0]  s_data_in_x;"<<endl;
    os<<"   output logic signed ["<<T-1<<":0]  m_data_out_y;"<<endl;
    os<<"   output logic               m_valid_y, s_ready_x;"<<endl;
    os<<"   logic                      wr_en_x, clear_acc, en_acc;"<<endl;
    os<<"   logic               ["<<addrbit_x-1<<":0]  addr_x;"<<endl;
    os<<"   logic               ["<<addrbit_f-1<<":0]  addr_f;"<<endl;
    os<<"   logic ["<<P-1<<"]["<<T-1<<"]               data_out_y"<<endl;
    os<<endl;
    os<<"   generate"<<endl;
    os<<"       genvar i;"<<endl;
    os<<"       for(i = 0; i < "<<P<<";i = i + 1) begin: control"<<endl;
    os<<"           Control control(s_valid_x, m_ready_y, s_ready_x, m_valid_y, clk, reset, wr_en_x, clear_acc, en_acc, addr_x, addr_f);"<<endl;
    os<<"       end"<<endl<<endl;
    os<<"   endgenerate"<<endl;
    os<<"   generate"<<endl;
    os<<"       genvar i;"<<endl;
    os<<"       for(i = 0; i < "<<P<<";i = i + 1) begin: control"<<endl;
    os<<"           Data_path data_path(m_data_out_y, addr_x, addr_f, wr_en_x, clear_acc, en_acc, clk, s_data_in_x, reset);"<<endl;
    os<<"       end"<<endl<<endl;
    os<<"   endgenerate"<<endl;
    os<<endl;
    os<<"endmodule"<<endl<<endl;


    os<<"module Control(s_valid_x, m_ready_y, s_ready_x, m_valid_y, clk, reset, wr_en_x, clear_acc, en_acc, addr_x, addr_f);"<<endl;
    os<<"   output logic ["<<addrbit_x-1<<":0] addr_x;"<<endl;
    os<<"   output logic ["<<addrbit_f-1<<":0] addr_f;"<<endl;
    os<<"   output logic       wr_en_x, clear_acc, en_acc;"<<endl;
    os<<"   input              s_valid_x, m_ready_y, clk, reset;"<<endl;
    os<<"   output logic       s_ready_x, m_valid_y;"<<endl;
    os<<"   logic              jishuqi;"<<endl;
    os<<"   logic              i; "<<endl;
    os<<"   logic              kk;"<<endl;
    os<<"   logic              kk_control;"<<endl;
    os<<"   logic              ii;"<<endl;
    os<<"   logic              jj;"<<endl;
    os<<"   logic              en_delay;"<<endl;
    os<<"   logic              en_delay2;"<<endl;
    os<<"   logic              count;"<<endl;
    os<<"   logic ["<<outbit<<":0]        count_output;"<<endl;
    os<<"   logic              validy_count;"<<endl;
    os<<"   logic              validy_en;"<<endl;
    os<<"   logic              validy_en2;"<<endl;
    os<<"   logic              validy_en3;"<<endl<<endl;
    os<<"   always_comb begin"<<endl;
    os<<"       if(s_valid_x == 1 && i == 0) begin"<<endl;
    os<<"           s_ready_x = 1;"<<endl;
    os<<"           wr_en_x = 1;"<<endl;
    os<<"       end"<<endl;
    os<<"       else begin"<<endl;
    os<<"           s_ready_x = 0;"<<endl;
    os<<"           wr_en_x = 0;"<<endl;
    os<<"       end"<<endl;
    os<<"   end"<<endl<<endl;
    os<<"   always_ff @(posedge clk) begin"<<endl;
    os<<"       if (addr_x == 0 && kk_control == 0) begin"<<endl;
    os<<"           kk <= kk + 1;"<<endl;
    os<<"           kk_control <= 1;"<<endl;
    os<<"       end"<<endl;
    os<<"	else if (addr_x != 0) begin"<<endl;
    os<<"           kk_control <= 0;"<<endl;
    os<<"	end"<<endl;
    os<<"       if (reset == 1) begin"<<endl;
    os<<"           addr_x <= 0;"<<endl;
    os<<"           addr_f <= "<<M-1<<";"<<endl;
    os<<"           m_valid_y <= 0;"<<endl;
    os<<"           jishuqi <= 0;"<<endl;
    os<<"           i <= 0;"<<endl;
    os<<"           kk <= 1; "<<endl;
    os<<"	    kk_control <= 0;"<<endl;
    os<<"           ii <= 0;"<<endl;
    os<<"           jj <= 0;"<<endl;
    os<<"           en_delay <= 0;"<<endl;
    os<<"           en_delay2 <= 0;"<<endl;
    os<<"           count <= 0;"<<endl;
    os<<"           count_output <= 0;"<<endl;
    os<<"           validy_count <= 0;"<<endl;
    os<<"       end"<<endl;
    os<<"       else begin"<<endl;
    os<<"           if(wr_en_x == 1 && addr_x < "<<N-1<<" && i == 0) begin"<<endl;
    os<<"               addr_x <=  addr_x + 1;"<<endl;
    os<<"               en_acc <= 0;"<<endl;
    os<<"               clear_acc <= 1;"<<endl;
    os<<"               count <= 0;"<<endl;
    os<<"               validy_en <= 0; "<<endl;
    os<<"               validy_en2 <= 0;"<<endl;
    os<<"               validy_en3 <= 0;"<<endl;
    os<<"           end"<<endl;
    os<<"       end"<<endl;
    os<<"       if (addr_x == "<<N-1<<" && wr_en_x == 1) begin"<<endl;
    os<<"           i <= 1;"<<endl;
    os<<"       end"<<endl;
    os<<"       if (i == 1)begin"<<endl;
    os<<"           if (jj == 0 && addr_x == "<<N-1<<" && addr_f == "<<M-1<<") begin"<<endl;
    os<<"               jishuqi <= 0;"<<endl;
    os<<"               addr_x <= 0;"<<endl;
    os<<"               addr_f <= 0;"<<endl;
    os<<"               clear_acc <= 0;"<<endl;
    os<<"               jj <= 1;"<<endl;
    os<<"           end"<<endl;
    os<<"           else if (addr_f <= "<<M-2<<") begin"<<endl;
    os<<"               en_acc <= 1;"<<endl;
    os<<"               jishuqi <= 0;"<<endl;
    os<<"               addr_x <= addr_x + 1;"<<endl;
    os<<"               addr_f <= addr_f + 1;"<<endl;
    os<<"           end"<<endl;
    os<<"       end"<<endl;
    os<<"       if (addr_f == "<<M-1<<" && kk >= 1 && count == 0 && i == 1) begin"<<endl;
    os<<"           en_delay <= 1;"<<endl;
    os<<"           count <= 1;"<<endl;
    os<<"       end"<<endl;
    os<<"       if (count == 1 && kk >= 1 && en_delay == 1) begin"<<endl;
    os<<"           en_delay2 <= 1;"<<endl;
    os<<"           en_delay <= 0;"<<endl;
    os<<"       end"<<endl;
    os<<"       if(en_delay2 == 1) begin"<<endl;
    os<<"           en_acc <= 0;"<<endl;
    os<<"           en_delay2 <= 0;"<<endl;
    os<<"       end"<<endl;
    os<<"       if (m_valid_y == 1 && m_ready_y == 1 && validy_count == 0) begin"<<endl;
    os<<"           count_output <= count_output + 1;"<<endl;
    os<<"           validy_count <= 1;"<<endl;
    os<<"       end"<<endl;
    os<<"       if(count_output == "<<N-M+1<<") begin"<<endl;
    os<<"           m_valid_y <= 0;"<<endl;
    os<<"           validy_count <= 0;"<<endl;
    os<<"           i <= 0;"<<endl;
    os<<"           count_output <= 0;"<<endl;
    os<<"           jj <= 0;"<<endl;
    os<<"           addr_x <= 0;"<<endl;
    os<<"           addr_f <= "<<M-1<<";"<<endl;
    os<<"           kk <= 1;"<<endl;
    os<<"       end"<<endl;
    os<<"       if (addr_x == "<<N-1<<" && addr_f == "<<M-1<<" && count_output == "<<N-M<<" && m_valid_y == 1 && m_ready_y == 1) begin"<<endl;
    os<<"           m_valid_y <= 0;"<<endl;
    os<<"       end"<<endl;
    os<<"       if (kk >= 1 && ii == 0 && i == 1 && addr_f == "<<M-1<<" && m_valid_y == 0) begin"<<endl;
    os<<"           validy_en <= 1;"<<endl;
    os<<"           ii <= 1;"<<endl;
    os<<"       end"<<endl;
    os<<"       if (validy_en == 1) begin"<<endl;
    os<<"           validy_en <= 0;"<<endl;
    os<<"           validy_en2 <= 1;"<<endl;
    os<<"       end"<<endl;
    os<<"       if (validy_en2 == 1) begin"<<endl;
    os<<"           m_valid_y <= 1;"<<endl;
    os<<"           validy_en2 <= 0;"<<endl;
    os<<"       end"<<endl;
    os<<"       if (count_output < "<<N-M<<" && m_valid_y == 1 && m_ready_y == 1 && i == 1) begin"<<endl;
    os<<"           m_valid_y <=  0;"<<endl;
    os<<"           validy_count <= 0;"<<endl;
    os<<"           clear_acc <= 1;"<<endl;
    os<<"           count <= 0;"<<endl;
    os<<"           addr_x <= addr_x - "<<M-2<<";"<<endl;
    os<<"           addr_f <= 0;"<<endl;
    os<<"       end"<<endl;
    os<<"       if (clear_acc == 1 && i == 1) begin"<<endl;
    os<<"           clear_acc <= 0;"<<endl;
    os<<"           ii <= 0;"<<endl;
    os<<"       end"<<endl;
    os<<"   end"<<endl<<endl;
    os<<"endmodule"<<endl<<endl;


    os<<"module Data_path(data_out, addr_x, addr_f, wr_en_x, clear_acc, en_acc, clk, s_data_in_x, reset);"<<endl;
   	long int i = 1;
	long int zhuangtai2 = -(i<<(T-1));
	long int zhuangtai3 = (i<<(T-1))-1;

    os<<"   output logic signed ["<<T-1<<":0] data_out;"<<endl;
    os<<"   logic signed ["<<T-1<<":0] data_out_Relu;"<<endl;
    os<<"   input        ["<<addrbit_x-1<<":0]        addr_x;"<<endl;
    os<<"   input        ["<<addrbit_f-1<<":0]        addr_f;"<<endl;
    os<<"   input                     wr_en_x, clear_acc, en_acc, clk, reset;"<<endl;
    os<<"   input signed ["<<T-1<<":0]       s_data_in_x;"<<endl;
    os<<"   logic signed ["<<T-1<<":0]       x_out;"<<endl;
    os<<"   logic signed ["<<T-1<<":0]       f_out;"<<endl;
    os<<"   logic signed ["<<2*T-1<<":0]       cheng;"<<endl;
    os<<"   logic signed ["<<2*T-1<<":0]       jia;"<<endl;
    os<<"   logic signed ["<<T-1<<":0]       pipe;"<<endl;
    os<<"   logic signed ["<<T-1<<":0]	     zhuangtai2;"<<endl;
    os<<"   logic signed ["<<T-1<<":0]	     zhuangtai3;"<<endl;
    os<<"   memory #("<<T<<", "<<N<<", "<<addrbit_x<<") mem_x(.clk(clk), .data_in(s_data_in_x), .data_out(x_out), .addr(addr_x), .wr_en(wr_en_x));"<<endl;
    os<<"   "<<romModName<<" rom"<<"(clk, addr_f, f_out);"<<endl<<endl;
    os<<"   always_comb begin"<<endl<<endl;
    os<<"	if (data_out_Relu <= 0) begin"<<endl;
    os<<"               data_out = 0"<<";"<<endl;
    os<<"	end"<<endl;
    os<<"	else if (data_out_Relu > 0) begin"<<endl;
    os<<"               data_out = data_out_Relu"<<";"<<endl;
    os<<"	end"<<endl<<endl;
    os<<"       if (en_acc == 1) begin"<<endl;
    os<<"           cheng = f_out * x_out;"<<endl;
    os<<"           jia = pipe + data_out_Relu;"<<endl;
    os<<"       end"<<endl;
    os<<"       else begin"<<endl;
    os<<"           cheng = 0;"<<endl;
    os<<"           jia = 0;"<<endl;
    os<<"       end"<<endl;
    os<<"   end"<<endl<<endl;
    os<<"   always_ff @(posedge clk) begin"<<endl;
    os<<"       if (clear_acc == 1 || reset == 1) begin"<<endl;
    os<<"           data_out_Relu <= 0;"<<endl;
    os<<"           pipe <= 0;"<<endl;
    os<<"	          zhuangtai2 <= "<<zhuangtai2<<";"<<endl;
    os<<"	          zhuangtai3 <= "<<zhuangtai3<<";"<<endl;
    os<<"       end"<<endl;
    os<<"       if (en_acc == 1) begin"<<endl;

    os<<"           if (cheng < zhuangtai2) begin"<<endl;
    os<<"               pipe <= zhuangtai2;"<<endl;
    os<<"           end"<<endl;
    os<<"           else if (cheng > zhuangtai3)begin"<<endl;
    os<<"               pipe <= zhuangtai3;"<<endl;
    os<<"           end"<<endl;
    os<<"           else begin"<<endl;
    os<<"               pipe <= cheng"<<";"<<endl;
    os<<"           end"<<endl<<endl<<endl;

    os<<"           if (jia < zhuangtai2) begin"<<endl;
    os<<"               data_out_Relu <= zhuangtai2;"<<endl;
    os<<"           end"<<endl;
    os<<"           else if (jia > zhuangtai3)begin"<<endl;
    os<<"               data_out_Relu <= zhuangtai3;"<<endl;
    os<<"           end"<<endl;
    os<<"           else begin"<<endl;
    os<<"               data_out_Relu <= jia"<<";"<<endl;
    os<<"           end"<<endl;

    os<<"       end"<<endl;
    os<<"   end"<<endl<<endl;
    os<<"endmodule"<<endl<<endl<<endl;


    os<<"module memory(clk, data_in, data_out, addr, wr_en);"<<endl;
    os<<"   parameter                          WIDTH = 8, SIZE = 256, LOGSIZE = 8;"<<endl;
    os<<"   input signed [WIDTH-1:0]           data_in;"<<endl;
    os<<"   output logic signed [WIDTH-1:0]    data_out;"<<endl;
    os<<"   input [LOGSIZE-1:0]                addr;"<<endl;
    os<<"   input                              clk, wr_en;"<<endl;
    os<<"   logic signed [SIZE-1:0][WIDTH-1:0] mem;	"<<endl<<endl;
    os<<"   always_ff @(posedge clk) begin"<<endl;
    os<<"       data_out <= mem[addr];"<<endl;
    os<<"       if(wr_en)"<<endl;
    os<<"           mem[addr] <= data_in;"<<endl;
    os<<"       end"<<endl;
    os<<"endmodule"<<endl<<endl;


    if (M > constVector.size()) {
        cout << "ERROR: constVector does not contain enough data for the requested design" << endl;
        cout << "The design parameters requested require " << M << " numbers, but the provided data only have " << constVector.size() << " constants" << endl;
        assert(false);
    }

    // Generate a ROM for f with constants in constVector, T bits, and the given name
    genROM(constVector, T, romModName, os);

}


// Part 3: Generate a hardware system with three layers interconnected.
// Layer 1: Input length: N, filter length: M1, output length: L1 = N-M1+1
// Layer 2: Input length: L1, filter length: M2, output length: L2 = L1-M2+1
// Layer 3: Input length: M2, filter length: M3, output length: L3 = L2-M3+1
// T is the number of bits
// A is the number of multipliers your overall design may use.
// Your goal is to build the highest-throughput design that uses A or fewer multipliers
// constVector holds all the constants for your system (all three layers, in order)
void genAllLayers(int N, int M1, int M2, int M3, int T, int A, vector<int>& constVector, string modName, ofstream &os) {

    // Here you will write code to figure out the best values to use for P1, P2, and P3, given
    // mult_budget.
    int P1 = 1; // replace this with your optimized value
    int P2 = 1; // replace this with your optimized value
    int P3 = 1; // replace this with your optimized value

    // output top-level module


    os << "module " << modName << "();" << endl;
    os << "   // this module should instantiate three convolution modules and wire them together" << endl;
    os << "endmodule" << endl;





    // -------------------------------------------------------------------------
    // Split up constVector for the three layers
    int start = 0;
    int stop = M1;
    vector<int> constVector1(&constVector[start], &constVector[stop]);

    // layer 2's W matrix is M2 x M1 and its B vector has size M2
    start = stop;
    stop = start+M2;
    vector<int> constVector2(&constVector[start], &constVector[stop]);

    // layer 3's W matrix is M3 x M2 and its B vector has size M3
    start = stop;
    stop = start+M3;
    vector<int> constVector3(&constVector[start], &constVector[stop]);

    if (stop > constVector.size()) {
        cout << "ERROR: constVector does not contain enough data for the requested design" << endl;
        cout << "The design parameters requested require " << stop << " numbers, but the provided data only have " << constVector.size() << " constants" << endl;
        assert(false);
    }
    // --------------------------------------------------------------------------


    // generate the three layer modules
    string subModName1 = "layer1_" + to_string(N) + "_" + to_string(M1) + "_" + to_string(T) + "_" + to_string(P1);
    genLayer(N, M1, T, P1, constVector1, subModName1, os);
    int L1 = N-M1+1;

    string subModName2 = "layer2_" + to_string(L1) + "_" + to_string(M2) + "_" + to_string(T) + "_" + to_string(P2);
    genLayer(L1, M2, T, P2, constVector2, subModName2, os);

    int L2 = L1-M2+1;
    string subModName3 = "layer3_" + to_string(L2) + "_" + to_string(M2) + "_" + to_string(T) + "_" + to_string(P3);
    genLayer(L2, M3, T, P3, constVector3, subModName3, os);

    // You will need to add code in the module at the top of this function to stitch together insantiations of these three modules

}



void printUsage() {
    cout << "Usage: ./gen MODE ARGS" << endl << endl;

    cout << "   Mode 1: Produce one convolution module (Part 1 and Part 2)" << endl;
    cout << "      ./gen 1 N M T P const_file" << endl;
    cout << "      See project description for explanation of parameters." << endl;
    cout << "      Example: produce a convolution with input vector of length 16, filter of length 4, parallelism 1" << endl;
    cout << "               and 16 bit words, with constants stored in file const.txt" << endl;
    cout << "                   ./gen 1 16 4 16 1 const.txt" << endl << endl;

    cout << "   Mode 2: Produce a system with three interconnected convolution module (Part 3)" << endl;
    cout << "      Arguments: N, M1, M2, M3, T, A, const_file" << endl;
    cout << "      See project description for explanation of parameters." << endl;
    cout << "              e.g.: ./gen 2 16 4 5 6 15 16 const.txt" << endl << endl;
}
