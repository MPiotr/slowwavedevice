#include<QtXml/QDomDocument>
#include<vector>

#ifndef vec
#define vec vector<double>
#endif
using namespace std;
bool setXMLEntry(QDomNode* node, char *name, double *par);
bool setXMLEntry(QDomNode* node, char *name, int *par);
bool setXMLEntry(QDomNode* node, char *name, char *par);
bool setXMLEntry(QDomDocument* doc, char *name, char *par);
bool setXMLEntry(QDomDocument* doc, char *name, double *par, vector<vec> *iterPar = NULL, vector<string> *iterNames = NULL);
bool setXMLEntry(QDomDocument* doc, char *name, int *par);
bool setXMLattribute(QDomDocument* doc, char *entryname, char* atributename, int *par);
bool setXMLattribute(QDomNode*  node, char* atributename, int *par);