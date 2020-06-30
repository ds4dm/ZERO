
.. _program_listing_file_src_support_utils.cpp:

Program Listing for File utils.cpp
==================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_support_utils.cpp>` (``src/support/utils.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #include "support/utils.h"
   #include <armadillo>
   #include <boost/log/trivial.hpp>
   #include <fstream>
   
   using namespace std;
   using namespace arma;
   
   arma::sp_mat Utils::resizePatch(const arma::sp_mat &mat, const unsigned int nR,
                                   const unsigned int nC) {
     arma::sp_mat mMat(nR, nC);
     mMat.zeros();
     if (nR >= mat.n_rows && nC >= mat.n_cols) {
        if (mat.n_rows >= 1 && mat.n_cols >= 1)
           mMat.submat(0, 0, mat.n_rows - 1, mat.n_cols - 1) = mat;
     } else {
        if (nR <= mat.n_rows && nC <= mat.n_cols)
           mMat = mat.submat(0, 0, nR, nC);
        else
           throw ZEROException(ZEROErrorCode::OutOfRange,
                               "Either both dimension should be smaller or larger.");
     }
     return mMat;
   }
   
   // For arma::mat
   arma::mat Utils::resizePatch(const arma::mat &mat, const unsigned int nR, const unsigned int nC) {
     arma::mat mMat(nR, nC);
     mMat.zeros();
     if (nR >= mat.n_rows && nC >= mat.n_cols) {
        if (mat.n_rows >= 1 && mat.n_cols >= 1)
           mMat.submat(0, 0, mat.n_rows - 1, mat.n_cols - 1) = mat;
     } else {
        if (nR <= mat.n_rows && nC <= mat.n_cols)
           mMat = mat.submat(0, 0, nR, nC);
        else
           throw ZEROException(ZEROErrorCode::OutOfRange,
                               "Either both dimension should be smaller or larger.");
     }
     return mMat;
   }
   
   // For arma::vec
   arma::vec Utils::resizePatch(const arma::vec &mat, const unsigned int nR) {
     arma::vec mMat(nR);
     mMat.zeros();
     if (nR > mat.n_rows)
        mMat.subvec(0, mat.n_rows - 1) = mat;
     else
        mMat = mat.subvec(0, nR);
     return mMat;
   }
   
   void Utils::appendSave(const sp_mat &matrix, 
                          const string out,     
                          const string header,  
                          bool erase            
                          )
   {
     // Using C++ file operations to copy the data into the target given by @out
     unsigned int nR{0}, nC{0}, nnz{0};
   
     ofstream outfile(out, erase ? ios::out : ios::app);
   
     nR  = matrix.n_rows;
     nC  = matrix.n_cols;
     nnz = matrix.n_nonzero;
   
     outfile << header << "\n";
     outfile << nR << "\t" << nC << "\t" << nnz << "\n";
     for (auto it = matrix.begin(); it != matrix.end(); ++it)
        outfile << it.row() << "\t" << it.col() << "\t" << (*it)
                << "\n"; // Write the required information of sp_mat
     outfile << "\n";
     outfile.close(); // and close it
   }
   
   long int Utils::appendRead(sp_mat &matrix,  
                              const string in, 
                              long int pos, 
                              const string header 
                              )
   {
     unsigned int nR = 0, nC = 0, nnz = 0;
   
     ifstream infile(in, ios::in);
     infile.seekg(pos);
   
     string headerCheckwith;
     infile >> headerCheckwith;
   
     if (header != "" && header != headerCheckwith)
        throw ZEROException(ZEROErrorCode::InvalidData,
                            "Wrong header. Expected " + header + " found " + headerCheckwith);
   
     infile >> nR >> nC >> nnz;
     if (nR == 0 || nC == 0)
        matrix.set_size(nR, nC);
     else {
        arma::umat locations(2, nnz);
        arma::vec values(nnz);
   
        unsigned int r = 0, c = 0;
        double val = 0;
   
        for (unsigned int i = 0; i < nnz; ++i) {
           infile >> r >> c >> val;
           locations(0, i) = r;
           locations(1, i) = c;
           values(i)       = val;
        }
        matrix = arma::sp_mat(locations, values, nR, nC);
     }
   
     pos = infile.tellg();
     infile.close();
   
     return pos;
   }
   
   void appendSave(const vector<double> v, const string out, const string header, bool erase) {
     ofstream outfile(out, erase ? ios::out : ios::app);
     outfile << header << "\n" << v.size() << "\n";
     for (const double x : v)
        outfile << x << "\n";
     outfile.close();
   }
   
   long int appendRead(vector<double> &v, const string in, long int pos, const string header) {
     unsigned long int size = 0;
     ifstream infile(in, ios::in);
     infile.seekg(pos);
     string headerCheckwith;
     infile >> headerCheckwith;
   
     if (header != "" && header != headerCheckwith)
        throw ZEROException(ZEROErrorCode::InvalidData,
                            "Wrong header. Expected " + header + " found " + headerCheckwith);
   
     infile >> size;
   
     v.resize(size);
     for (unsigned int i = 0; i < size; ++i)
        infile >> v[i];
     pos = infile.tellg();
     infile.close();
     return pos;
   }
   
   void Utils::appendSave(const vec &matrix,   
                          const string out,    
                          const string header, 
                          bool erase           
   ) {
     // Using C++ file operations to copy the data into the target given by @out
     unsigned int nR{0};
   
     ofstream outfile(out, erase ? ios::out : ios::app);
   
     nR = matrix.n_rows;
   
     outfile << header << "\n";
   
     outfile << nR << "\n";
     for (auto it = matrix.begin(); it != matrix.end(); ++it)
        outfile << (*it) << "\n"; // Write the required information of sp_mat
     outfile << "\n";
     outfile.close(); // and close it
   }
   
   long int Utils::appendRead(vec &matrix,     
                              const string in, 
                              long int pos, 
                              const string header 
   ) {
     unsigned int nR;
     string buffers;
     string checkwith;
     ifstream inFile(in, ios::in);
     inFile.seekg(pos);
   
     inFile >> checkwith;
     if (header != "" && checkwith != header)
        throw ZEROException(ZEROErrorCode::InvalidData,
                            "Wrong header. Expected " + header + " found " + checkwith);
     inFile >> nR;
     matrix.zeros(nR);
     for (unsigned int i = 0; i < nR; ++i) {
        double val;
        inFile >> val;
        matrix.at(i) = val;
     }
   
     pos = inFile.tellg();
     inFile.close();
   
     return pos;
   }
   
   void Utils::appendSave(const long int v, const string out, const string header, bool erase)
   {
     ofstream outfile(out, erase ? ios::out : ios::app);
     outfile << header << "\n";
     outfile << v << "\n";
     outfile.close();
   }
   
   long int Utils::appendRead(long int &v, const string in, long int pos, const string header) {
     ifstream infile(in, ios::in);
     infile.seekg(pos);
   
     string headerCheckwith;
     infile >> headerCheckwith;
   
     if (header != "" && header != headerCheckwith)
        throw ZEROException(ZEROErrorCode::InvalidData,
                            "Wrong header. Expected " + header + " found " + headerCheckwith);
   
     long int val;
     infile >> val;
     v = val;
   
     pos = infile.tellg();
     infile.close();
   
     return pos;
   }
   
   void Utils::appendSave(const unsigned int v, const string out, const string header, bool erase)
   {
     ofstream outfile(out, erase ? ios::out : ios::app);
     outfile << header << "\n";
     outfile << v << "\n";
     outfile.close();
   }
   
   long int Utils::appendRead(unsigned int &v, const string in, long int pos, const string header) {
     ifstream infile(in, ios::in);
     infile.seekg(pos);
   
     string headerCheckwith;
     infile >> headerCheckwith;
   
     if (header != "" && header != headerCheckwith)
        throw ZEROException(ZEROErrorCode::InvalidData,
                            "Wrong header. Expected " + header + " found " + headerCheckwith);
   
     unsigned int val;
     infile >> val;
     v = val;
   
     pos = infile.tellg();
     infile.close();
   
     return pos;
   }
   
   void Utils::appendSave(const string v, const string out, bool erase)
   {
     ofstream outfile(out, erase ? ios::out : ios::app);
     outfile << v << "\n";
     outfile.close();
   }
   
   long int Utils::appendRead(string &v, const string in, long int pos) {
     ifstream infile(in, ios::in);
     infile.seekg(pos);
   
     string val;
     infile >> val;
     v = val;
   
     pos = infile.tellg();
     infile.close();
   
     return pos;
   }
   unsigned long int Utils::vecToNum(std::vector<short int> binary) {
     unsigned long int number = 0;
     unsigned int posn        = 1;
     while (!binary.empty()) {
        short int bit = (binary.back() + 1) / 2; // The least significant bit
        number += (bit * posn);
        posn *= 2;         // Update place value
        binary.pop_back(); // Remove that bit
     }
     return number;
   }
   
   std::vector<short int> Utils::numToVec(unsigned long int number, const unsigned long nCompl) {
     std::vector<short int> binary{};
     for (unsigned int vv = 0; vv < nCompl; vv++) {
        binary.push_back(number % 2);
        number /= 2;
     }
     std::for_each(binary.begin(), binary.end(), [](short int &vv) { vv = (vv == 0 ? -1 : 1); });
     std::reverse(binary.begin(), binary.end());
     return binary;
   }
   
   bool Utils::containsConstraint(arma::sp_mat &A, const vec &b, const arma::vec &lhs,
                                  const double &rhs, const double tol) {
     if (lhs.size() != A.n_cols)
        return false;
     for (int i = 0; i < A.n_rows; ++i) {
        bool res = true;
        for (int j = 0; j < A.n_cols; ++j) {
           if (std::abs(lhs.at(j) - A.at(i, j)) > tol) {
             res = false;
             break;
           }
        }
        if (res && std::abs(b.at(i) - rhs) < tol) {
           return true;
        }
     }
     return false;
   }
   
   bool Utils::containsElement(const vec &b, const double &element, const double tol) {
     for (unsigned int i = 0; i < b.size(); ++i) {
        if (std::abs(b.at(i) - element) < tol)
           return true;
     }
     return false;
   }
   
   bool Utils::containsRow(const sp_mat &A, const arma::vec &row, const double tol) {
   
     if (row.size() != A.n_cols)
        return false;
     for (int i = 0; i < A.n_rows; ++i) {
        bool res = true;
        for (int j = 0; j < A.n_cols; ++j) {
           if (std::abs(row.at(j) - A.at(i, j)) > tol) {
             res = false;
             break;
           }
        }
        if (res)
           return true;
     }
     return false;
   }
   bool Utils::containsConstraint(sp_mat &A, const vec &b, const sp_mat &lhs, const double &rhs,
                                  const double tol) {
     if (lhs.n_rows > 1)
        return false;
     arma::vec Ai = arma::vec{lhs};
     return Utils::containsConstraint(A, b, Ai, rhs, tol);
   }
   arma::vec Utils::normalize(const arma::vec input) {
     double max  = input.max();
     double min  = std::abs(input.min());
     double norm = 0;
   
     if (max > min)
        norm = max;
     else
        norm = min;
   
     return input / norm;
   }
