//   Copyright (c)  2010  Anna Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include "OnlineHelp.H"
#include "CoCoA/error.H"
//#include "CoCoA/io.H"

#include <algorithm>
using std::transform;
#include <cctype>
//using std::tolower;
#include <fstream>
using std::ifstream;
#include <iostream>
using std::clog;
using std::endl;
#include <string>
using std::string;
#include <vector>
using std::vector;


extern string packageDir;

namespace CoCoA
{
namespace OnlineHelp
{

  class index
  {
  public: 
    class entry
    {
    public:
      entry(const std::string& title): myTitleValue(title) {}
      // ~entry(); // default is OK
      const std::string& myTitle() const {return myTitleValue;}
      const std::vector<std::string>& myKeys() const {return myKeysValue;}
      const std::vector<std::string>& myTypes() const {return myTypesValue;}
      const std::vector<std::string>& myRtnTypes() const {return myRtnTypesValue;}
      const std::vector<std::string>& mySeeAlso() const {return mySeeAlsoValue;}
      const std::string& mySyntax() const {return mySyntaxValue;}
      void myAddKey(const std::string& key);
      void myAddType(const std::string& key);
      void myAddRtnType(const std::string& key);
      void myAddSee(const std::string& key);
      void myAddKeys(std::ifstream& in);
      void myAddTypes(std::ifstream& in);
      void myAddSeeAlso(std::ifstream& in);
      void myAddSyntax(std::ifstream& in);

    private: // data members
      std::string myTitleValue;
      std::vector<std::string> myKeysValue;
      std::vector<std::string> myTypesValue;
      std::vector<std::string> myRtnTypesValue;
      std::vector<std::string> mySeeAlsoValue;
      std::string mySyntaxValue;
    };
 
  public:
    index();
    // ~index(); // default is OK
    void myLoad(std::ostream &out);
    const std::string& myTitle(std::size_t i) const {return myVec[i].myTitle();}
    const std::vector<std::string>& myKeys(std::size_t i) const {return myVec[i].myKeys();}
    const std::vector<std::string>& myTypes(std::size_t i) const {return myVec[i].myTypes();}
    const std::vector<std::string>& myRtnTypes(std::size_t i) const {return myVec[i].myRtnTypes();}
    const std::vector<std::string>& mySeeAlso(std::size_t i) const {return myVec[i].mySeeAlso();}
    const std::string& mySyntax(std::size_t i) const {return myVec[i].mySyntax();}
    std::size_t mySize() const {return myVec.size();}
    std::string myManCoCoAVersion() const {return myManCoCoAVersionValue;}
    std::string myManDate() const {return myManDateValue;}

  private:
    void myAddEntry(std::ifstream& in, const std::string& title, const std::string& ClosingTag);
  
  private: // data members
    std::vector<entry> myVec;
    std::string myManDateValue;
    std::string myManCoCoAVersionValue;
  };

  // ----------------------------------------------------------------------
  // fwd decl -- defined later in this file
  std::string LowerCase(const std::string& str);
  bool IsStr1BeforeStr2(std::ifstream& in, const std::string& s1, const std::string& s2);
  std::string StringInTag(const std::string& line, const std::string& XMLTag);
  std::string SkipToStringInTag(std::ifstream& in, const std::string& XMLTag);
  void PrintCommAndFunc(std::string& line, std::ostream &out);
  void PrintCommAndFuncRtn(std::string& line, std::ostream &out);
  void PrintSeeAlso(const std::string& s, std::ostream &out);
  void PrintSyntax(const std::string& s, std::ostream &out);
  std::string CleanupLine(const std::string& line);
  std::string CleanupKeyword(const std::string& key);
  void PrintAllMatchesFor(const std::string& str, std::ostream &out);
  void PrintAllMatchesFor(const std::string& str,
                          const std::vector<string>& found,
                          std::ostream &out);

  // ----------------------------------------------------------------------
  inline bool IsSubstring(const std::string& s, const std::string& line) 
  { return line.find(s) != string::npos; }

  // ----------------------------------------------------------------------

  //-- index --

  //  const index& UniqueIndex();
  index& UniqueIndex();

  index::index()
  {
    myLoad(clog);
  }


  void index::myLoad(std::ostream &out)
  {
    out << "Loading CoCoAManual index ..." << std::flush;
    ifstream in(XMLFileName());
    if (in == 0)
      CoCoA_ERROR("Cannot find input file `CoCoAHelp.xml`", "OnlineHelp::index()");
    string line; // used only inside while loop
    string LastCommand;
    string NewCommand = "a";
    string ManCoCoAVersion;
    string ManDate;
    myVec.clear();
    while (!in.eof())
    {
      getline(in, line);
      if (IsSubstring("<!--", line)) continue;
      if (IsSubstring("<version>", line))
      {
        myManCoCoAVersionValue = SkipToStringInTag(in, "<cocoa_version>");
        myManDateValue = SkipToStringInTag(in, "<date>");
      }
      if (IsSubstring("<command>", line))
      {
        LastCommand = NewCommand;
        NewCommand = SkipToStringInTag(in, "<title>");
        // check for proper sorting
        if (LowerCase(LastCommand).compare(LowerCase(NewCommand))>0)
          std::cout << "OnlineHelp: unsorted entries: "
                    << LastCommand << " -- " << NewCommand << std::endl;
        myAddEntry(in, CleanupLine(NewCommand), "</command>");
      }
      if (IsSubstring("<section>", line))
        myAddEntry(in, CleanupLine(SkipToStringInTag(in, "<title>")), "</section>");
    }
    out << "... CoCoAManual index loaded" << endl;
  }


  void index::myAddEntry(std::ifstream& in, const std::string& title, const std::string& ClosingTag)
  {
    myVec.push_back(entry(title));
    entry& e(myVec.back());
    e.myAddKey(title);
    string line;
    getline(in, line);
    while (!in.eof() && !IsSubstring(ClosingTag, line))
    {
      if (IsSubstring("<keys>", line)) e.myAddKeys(in);
      if (IsSubstring("<types>", line)) e.myAddTypes(in);
      if (IsSubstring("<seealso>", line)) e.myAddSeeAlso(in);
      if (IsSubstring("<syntax>", line)) e.myAddSyntax(in);
      getline(in, line);
    }
  }


  //-- index::entry --
  void index::entry::myAddKey(const std::string& key)
  { myKeysValue.push_back(CleanupLine(LowerCase(key))); }

  void index::entry::myAddType(const std::string& t)
  { myTypesValue.push_back(t); }

  void index::entry::myAddRtnType(const std::string& t)
  { myRtnTypesValue.push_back(t); }

  void index::entry::myAddSee(const std::string& t)
  { mySeeAlsoValue.push_back(t); }

  void index::entry::myAddKeys(std::ifstream& in)
  {
    string line;
    getline(in, line);
    string s; // used only inside while loop
    while ((!in.eof()) && (!IsSubstring("</keys>",line)))
    {
      if ((s=StringInTag(line, "<key>")) != "") myAddKey(s);
      getline(in, line);
    }
  }

  void index::entry::myAddTypes(std::ifstream& in)
  {
    string line;
    getline(in, line);
    string s; // used only inside while loop
    while ((!in.eof()) && (!IsSubstring("</types>",line)))
    {
      if ((s=StringInTag(line, "<type>")) != "") myAddType(s);
      getline(in, line);
    }
  }


  void index::entry::myAddSeeAlso(std::ifstream& in)
  {
    string line;
    getline(in, line);
    string s; // used only inside while loop
    while ((!in.eof()) && (!IsSubstring("</seealso>",line)))
    {
      if ((s=StringInTag(line, "<see>")) != "") myAddSee(s);
      getline(in, line);
    }
  }


  void index::entry::myAddSyntax(std::ifstream& in)
  {
    string line;
    getline(in, line);
    while ((!in.eof()) && (!IsSubstring("</syntax>",line)))
    {
      if (IsSubstring("<type>", line))
        myAddType(StringInTag(line, "<type>"));
      if (IsSubstring("<rtn>", line))
        myAddRtnType(StringInTag(line, "<rtn>"));
      mySyntaxValue += "--> ";
      mySyntaxValue += CleanupLine(line);
      mySyntaxValue += "\n";
      getline(in, line);
    }
  }


  //----------------------------------------------------------------

  std::string LowerCase(const std::string& str)
  {
    string s(str); // slightly wasteful, but readable
    transform(s.begin(), s.end(), s.begin(), tolower);
    return s;
  }


  bool IsStr1BeforeStr2(std::ifstream& in, const std::string& s1, const std::string& s2)
  {
    string line;
    getline(in, line);
    while (!in.eof())
    {
      if ( IsSubstring(s2, line) ) return false;
      if ( IsSubstring(s1, line) ) return true;
      getline(in, line);
    }
    CoCoA_ERROR(CoCoA::ERR::SERIOUS, "IsStr1BeforeStr2");
    return false;
  }
  

  // only for opening and closing tag in the same line
  std::string StringInTag(const std::string& line, const std::string& XMLTag)
  {
    size_t open;
    if ( (open=line.find(XMLTag)) == string::npos) return "";
    
    string ClosedXMLTag = XMLTag;
    ClosedXMLTag.replace(0, 1, "</");
    size_t close = line.find(ClosedXMLTag);
    if (close == string::npos)
      CoCoA_ERROR(XMLTag+" closing tag not found in this line", "StringInTag");
    size_t StartPos = open + XMLTag.length();
    return line.substr(StartPos, close-StartPos);
  }


  std::string SkipToStringInTag(std::ifstream& in, const std::string& XMLTag)
  {
    string line;
    getline(in, line);
    string s;
    while (!in.eof())
    {
      if ( (s=StringInTag(line, XMLTag)) != "") return s;
      getline(in, line);
    }
    return "";
  }
  

  void ReplaceWith(std::string& line, const std::string& SFrom, const std::string& STo)
  {
    size_t open;
    while ( (open=line.find(SFrom)) != string::npos)
      line.replace(open, SFrom.length(), STo);
  }
  

  void ReplaceWithQuotes(std::string& line, const std::string& s)
  { ReplaceWith(line, s, "\""); }
  

  std::string& ProcessLine(std::string& line)
  {
    ReplaceWithQuotes(line, "<quotes>"); ReplaceWithQuotes(line, "</quotes>");
    ReplaceWithQuotes(line, "<tt>");     ReplaceWithQuotes(line, "</tt>");
    ReplaceWithQuotes(line, "<code>");   ReplaceWithQuotes(line, "</code>");
    ReplaceWithQuotes(line, "<ttref>");  ReplaceWithQuotes(line, "</ttref>");
    ReplaceWithQuotes(line, "<ref>");    ReplaceWithQuotes(line, "</ref>");
    ReplaceWith(line, "<em>", "*");      ReplaceWith(line, "</em>", "*");
    ReplaceWith(line, "<i>", "");        ReplaceWith(line, "</i>", "");
    ReplaceWith(line, "<verbatim>", ""); ReplaceWith(line, "</verbatim>", "");
    ReplaceWith(line, "<formula>", " "); ReplaceWith(line, "</formula>", " ");
    ReplaceWith(line, "&lt;", "<");
    ReplaceWith(line, "&gt;", ">");
    ReplaceWith(line, "&apos;", "'");
    ReplaceWith(line, "&amp;", "&");
    ReplaceWith(line, "<less_eq/>", "<=");
    ReplaceWith(line, "<par/>", "");
    ReplaceWith(line, "<par />", "");
    ReplaceWith(line, "<cocoa_version/>", UniqueIndex().myManCoCoAVersion());
    ReplaceWith(line, "<cocoa_date/>", UniqueIndex().myManDate());
    return line;
  }


  std::string CleanupLine(const std::string& line)
  {
    // CleanupLine is called by the ctor UniqueIndex
    string s(line);
    ReplaceWithQuotes(s, "<quotes>"); ReplaceWithQuotes(s, "</quotes>");
    ReplaceWithQuotes(s, "<tt>");     ReplaceWithQuotes(s, "</tt>");
    ReplaceWith(s, "<type>", "");     ReplaceWith(s, "</type>", "");
    ReplaceWith(s, "<rtn>", "");      ReplaceWith(s, "</rtn>", "");
    ReplaceWith(s, "<em>", "*");      ReplaceWith(s, "</em>", "*");
    ReplaceWith(s, "&lt;", "<");
    ReplaceWith(s, "&gt;", ">");
    ReplaceWith(s, "&apos;", "'");
    ReplaceWith(s, "&amp;", "&");
    ReplaceWith(s, "<less_eq/>", "<=");
    return s;
  }


  void SkipComment(std::ifstream& in, std::string& line)
  {
    while (  (!in.eof()) && !IsSubstring("-->", line) )  getline(in, line);
    if ( in.eof() )
      CoCoA_ERROR("missing comment close tag: eof", "SkipComment");
    const size_t close = line.find("-->");
    if (close != line.size()-3)
      CoCoA_ERROR("text after closed comment: \""+line+"\"", "SkipComment");
  }


  void PrintExample(std::ifstream& in, std::ostream &out)
  {
    out << "------<  example  >------" << endl;
    string line;
    getline(in, line);
    while ( !IsSubstring("</example>", line) )
    {
      out << CleanupLine(line) << endl;
      getline(in, line);
    }
    out << "------< end example >------" << endl;
  }
  

  void PrintExampleWithoutOutput(std::ifstream& in, std::ostream &out)
  {
    out << "------<  example  >------" << endl;
    string line;
    getline(in, line);
    while ( !IsSubstring("</example>", line) )
    {
      if (IsSubstring("/**/", line))  out << CleanupLine(line) << endl;
      getline(in, line);
    }
    out << "------< end example >------" << endl;
  }
  

  void PrintSyntax(const std::string& s, std::ostream &out)
  {
    const index& index(UniqueIndex());
    bool IsFirst=true;
    for (size_t i=0; i<index.mySize(); ++i)
      if ( LowerCase(index.myTitle(i)) == LowerCase(s) && !index.mySyntax(i).empty())
      {
        if (IsFirst)
        {
          IsFirst = false;
          //          out << "--====<  syntax  >====--" << endl;
        }
        out << index.mySyntax(i) << endl;
        break;
      }
  }
  

void PrintSearchExact(std::string title, std::ostream &out)
{
  ifstream in(XMLFileName());
  if (in == 0)
  {
    std::cout << "Cannot find input file `CoCoAHelp.xml`.  Aborting." << endl;
    abort();
  }
  string s = SkipToStringInTag(in, "<title>");
  while (!in.eof() && s != title)  s = SkipToStringInTag(in, "<title>");
  if (in.eof())
    CoCoA_ERROR("Not found exact title: " + title, "PrintSearchExact");
  out << "--============[ " << CleanupLine(title) << " ]=============--" << endl;
  PrintSyntax(title, out);
  string line;
  getline(in, line);
  while (!IsSubstring("<description>", line))  getline(in, line);
  //  out << "--====<  description  >====--" << endl;
  getline(in, line);
  while (!IsSubstring("</description>", line))
  {
    if ( IsSubstring("<obsolete_functions>", line))
      PrintAllMatchesFor("[obsolete]", out);
    else
    {
    if ( IsSubstring("<obsolescent_functions>", line))
      PrintAllMatchesFor("[obsolescent]", out);
    else
    {
    if ( IsSubstring("<commands_and_functions_for", line))
      PrintCommAndFunc(line, out);
    else
    {
      if ( IsSubstring("<commands_and_functions_rtn", line))
        PrintCommAndFuncRtn(line, out);
      else
      {
        if ( IsSubstring("<example>", line))  PrintExample(in, out);
        else
        {
          if ( IsSubstring("<!--", line))  SkipComment(in, line);
          else
            out << ProcessLine(line) << endl;
        }
      }
    }
    }
    }
    getline(in, line);
  }
  getline(in, line);
  while (!IsSubstring("<", line))   getline(in, line);
  PrintSeeAlso(title, out);
}


#if 0 
  // hopefully useless
  std::size_t SkipToLAngle(std::ifstream& in, std::string& line, std::size_t from)
  {
    size_t pos;
    while (!in.eof())
    {
      if ( (pos=line.find("<")) != string::npos) return pos;
      getline(in, line);
    }
    return string::npos;
  }
#endif


//----------------------------------------------------------------------

void SearchMatch(std::vector<std::string>& MatchingTitles, const std::string& s)
{
  const string ls(LowerCase(s));
  
  const index& index(UniqueIndex());
  for (size_t i=0; i<index.mySize(); ++i)
    for (size_t j=0; j<index.myKeys(i).size(); ++j)
      if (IsSubstring(ls, index.myKeys(i)[j]))
      {
        MatchingTitles.push_back(index.myTitle(i));
        break;
      }
}


void SearchMatch(std::vector<long>& MatchingEntries, const std::string& s)
{
  const string ls(LowerCase(s));
  
  const index& index(UniqueIndex());
  for (size_t i=0; i<index.mySize(); ++i)
    for (size_t j=0; j<index.myKeys(i).size(); ++j)
      if (IsSubstring(ls, index.myKeys(i)[j]))
      {
        MatchingEntries.push_back(i);
        break;
      }
}


void SearchType(std::vector<std::string>& MatchingTitles, const std::string& s)
{
  const index& index(UniqueIndex());
  for (size_t i=0; i<index.mySize(); ++i)
    for (size_t j=0; j<index.myTypes(i).size(); ++j)
      if (s == index.myTypes(i)[j])
      {
        MatchingTitles.push_back(index.myTitle(i));
        break;
      }
}


void SearchRtnType(std::vector<std::string>& MatchingTitles, const std::string& s)
{
  const index& index(UniqueIndex());
  for (size_t i=0; i<index.mySize(); ++i)
    for (size_t j=0; j<index.myRtnTypes(i).size(); ++j)
      if (s == index.myRtnTypes(i)[j])
      {
        MatchingTitles.push_back(index.myTitle(i));
        break;
      }
}


  void PrintCommAndFunc(std::string& TypeLine, std::ostream &out)
  {
    ReplaceWith(TypeLine, "<commands_and_functions_for type=\"", "");
    ReplaceWith(TypeLine, "\"></commands_and_functions_for>", "");
    std::vector<std::string> found;
    SearchType(found, TypeLine);
    for (size_t i=0; i<found.size(); ++i)
    {
      ifstream in(XMLFileName());  // if (in == 0) ...
      string s = SkipToStringInTag(in, "<title>");
      while (!in.eof() && s != found[i])  s = SkipToStringInTag(in, "<title>");
      if (in.eof())
        CoCoA_ERROR("Not found exact title: " + found[i], "PrintCommAndFunc");
      out << "? " << s << " -- ";
      string line;
      getline(in, line);
      while (!IsSubstring("<short_description>", line))  getline(in, line);
      out << StringInTag(line, "<short_description>") << endl;
    }
    out << "--==================<>===================--" << endl;
  }  


  void PrintCommAndFuncRtn(std::string& TypeLine, std::ostream &out)
  {
    ReplaceWith(TypeLine, "<commands_and_functions_rtn type=\"", "");
    ReplaceWith(TypeLine, "\"></commands_and_functions_rtn>", "");
    std::vector<std::string> found;
    SearchRtnType(found, TypeLine);
    for (size_t i=0; i<found.size(); ++i)
    {
      ifstream in(XMLFileName());  // if (in == 0) ...
      string s = SkipToStringInTag(in, "<title>");
      while (!in.eof() && s != found[i])  s = SkipToStringInTag(in, "<title>");
      if (in.eof())
        CoCoA_ERROR("Not found exact title: " + found[i], "PrintCommAndFuncRtn");
      out << "? " << s << " -- ";
      string line;
      getline(in, line);
      while (!IsSubstring("<short_description>", line))  getline(in, line);
      out << StringInTag(line, "<short_description>") << endl;
    }
    out << "--==================<>===================--" << endl;
  }  


  string CleanupKeyword(const std::string& key)
  {
    using std::isspace;
    string keyword(key);
    // remove trailing comments, spaces, ";"
    if (IsSubstring("--",keyword))  keyword.resize(keyword.find("--"));
    if (IsSubstring("//",keyword))  keyword.resize(keyword.find("//"));
    for (long i=keyword.length()-1; i>=0; --i)
      if (keyword[i]==' ' || keyword[i]=='\t' || keyword[i]==';')
        keyword.resize(i);
      else break;
    // Strip initial spaces
    const long len_keyword = keyword.size();
    long i=0;
    while (i < len_keyword && isspace(keyword[i])) ++i;
    keyword = keyword.substr(i); // discard initial whitespace
    return keyword;
  }
  

  void PrintMan(std::string keyword, std::ostream &out)
  {
    const index& index(UniqueIndex());
    keyword = CleanupKeyword(keyword);
    enum { SingleQuery, DoubleQuery } HelpType = SingleQuery;

    // Check whether there is a 2nd '?'; if so, skip it & whitespace
    if (keyword[0]=='?')
    {
      HelpType = DoubleQuery;
      keyword = CleanupKeyword(keyword.substr(1)); // skip fist char: '?'
    }

    vector<long> AllMatches;
    SearchMatch(AllMatches, keyword);
    if (AllMatches.empty())
    {
      out << "\nNo matches for \"" << keyword << "\"" << endl;
      return;
    }
    const long NumAllMatches = AllMatches.size();
    if (HelpType == SingleQuery && NumAllMatches == 1)
    {
      PrintSearchExact(index.myTitle(AllMatches[0]), out);
      return;
    }
    long ExactMatchPosn = -1;
    vector<long> ExactMatchKey;
    if (HelpType == SingleQuery)
      for (long i=0; i<NumAllMatches; ++i)
      {
        if (LowerCase(index.myTitle(AllMatches[i])) == LowerCase(keyword))
        {
          ExactMatchPosn = i;
          PrintSearchExact(index.myTitle(AllMatches[ExactMatchPosn]), out);
          break;
        }
        const std::vector<std::string>& keys_i = index.myKeys(AllMatches[i]);
        for (long n=keys_i.size()-1; n>=0; --n)
          if (LowerCase(keys_i[n]) == LowerCase(keyword))
          {
            ExactMatchKey.push_back(i);
            break;
          }
      }
    if (ExactMatchPosn == -1)
    {
      if (ExactMatchKey.size()==1)
      {
        ExactMatchPosn = ExactMatchKey[0];
        PrintSearchExact(index.myTitle(AllMatches[ExactMatchPosn]), out);
      }
      else
      {
        if (HelpType == SingleQuery)
          out << "--====<  No entry for \"" << keyword << "\" >====--";
        PrintAllMatchesFor(keyword, out);
        return;
      }
    }
    if (NumAllMatches > 6)
    {
      out << "\nTo see all " << NumAllMatches << " matches for \"" << keyword << "\":\n";
      out << " ?? " << keyword << endl;
    }
    else
    {
      if (ExactMatchPosn == -1)
      {
        out << "\nAll further matches for \"" << keyword << "\":\n";
        for (int i=0; i<NumAllMatches; ++i)
          out << " ? " << index.myTitle(AllMatches[i]) << endl;
      }
      else
      {
        vector<long> FurtherMatches;
        const vector<string>& EMSeeAlso = index.mySeeAlso(AllMatches[ExactMatchPosn]);
        for (int i=0; i<NumAllMatches; ++i)
        {
          bool AlredyCited = false;
          if (i == ExactMatchPosn) AlredyCited = true;
          for (size_t s=0; s<EMSeeAlso.size(); ++s)
            if ( index.myTitle(AllMatches[i]) == EMSeeAlso[s] )
              AlredyCited = true;
          if (!AlredyCited) FurtherMatches.push_back(AllMatches[i]);
        }
        if (!FurtherMatches.empty())
        {
          out << "\nAll " << FurtherMatches.size() << " further matches for \"" << keyword << "\":\n";
          for (size_t i=0; i<FurtherMatches.size(); ++i)
            out << " ? " << index.myTitle(FurtherMatches[i]) << endl;
        }
      }
    }
  }


void PrintAllMatchesFor(const std::string& str,
                        const std::vector<string>& found,
                        std::ostream &out)
  {
    if (found.empty())
    {
      out << "\nNo matches for \"" << str << "\"" << endl;
      return;
    }
    const int NumFound = found.size();
    out << "\nAll " << NumFound << " matches for \"" << str << "\":\n";
//    out << " (" << NumFound << ")" << endl;
    for (int i=0; i<NumFound; ++i)
      out << " ? " << found[i] << endl;
  }
  

  void PrintAllMatchesFor(const std::string& s, std::ostream &out)
  {
    vector<string> found;
    SearchMatch(found, s);
    PrintAllMatchesFor(s, found, out);
  }
  

  void ReloadMan(std::ostream &out)
  {
    UniqueIndex().myLoad(out);
  }
  

  void PrintSeeAlso(const std::string& s, std::ostream &out)
  {
    const index& index(UniqueIndex());
    for (size_t i=0; i<index.mySize(); ++i)
      if ( LowerCase(index.myTitle(i)) == LowerCase(s) )
      {
        if (index.mySeeAlso(i).size() != 0)
          out << "\nSee also:" << endl;
        for (size_t j=0; j<index.mySeeAlso(i).size(); ++j)
          out << " ? " << index.mySeeAlso(i)[j] << endl;
        break;
      }
  }
  

void PrintAllExamples(std::ostream &out)
{
  ifstream in(XMLFileName());
  if (in == 0)
  {
    std::cout << "Cannot find input file `CoCoAHelp.xml`.  Aborting." << endl;
    abort();
  }
  string s = "";
  while (!in.eof())
  {
    s = SkipToStringInTag(in, "<title>");
    if (in.eof()) break;
    out << "PrintLn \"-- " << s << " =============================\";" << endl;
    string line;
    getline(in, line);
    while (!IsSubstring("<description>", line))  getline(in, line);
    getline(in, line);
    while (!IsSubstring("</description>", line))
    {
      if ( IsSubstring("<example>", line) && !IsSubstring("<!--", line) )
        PrintExample(in, out);
      getline(in, line);
    }
  }
  out << "--- PrintAllExamples: done ---" << endl;
}


void PrintAllExamplesWithoutOutput(std::ostream &out)
{
  ifstream in(XMLFileName());
  if (in == 0)
  {
    std::cout << "Cannot find input file `CoCoAHelp.xml`.  Aborting." << endl;
    abort();
  }
  string s = "";
  while (!in.eof())
  {
    s = SkipToStringInTag(in, "<title>");
    if (in.eof()) break;
    out << "PrintLn \"-- " << s << " =============================\";" << endl;
    string line;
    getline(in, line);
    while (!IsSubstring("<description>", line))  getline(in, line);
    getline(in, line);
    while (!IsSubstring("</description>", line))
    {
      if ( IsSubstring("<example>", line) && !IsSubstring("<!--", line) )
        PrintExampleWithoutOutput(in, out);
      getline(in, line);
    }
  }
  out << "--- PrintAllExamples: done ---" << endl;
}


void PrintWordlist(std::ostream &out)
{
  ifstream in(XMLFileName());
  if (in == 0)
  {
    std::cout << "Cannot find input file `CoCoAHelp.xml`.  Aborting." << endl;
    abort();
  }
  string s = "";
  string line;
  getline(in, line);
  while (!in.eof())
  {
    while (!IsSubstring("<command>", line) && !in.eof())  getline(in, line);
    out << SkipToStringInTag(in, "<title>") << endl;
    getline(in, line);
  }
}


const char* XMLFileName()
{
  static const string UniqueCopy(packageDir+"/../CoCoAManual/CoCoAHelp.xml");
  return UniqueCopy.c_str();
}


//const index& UniqueIndex()
index& UniqueIndex()
{
  static index UniqueCopy;
  return UniqueCopy;
}

} // namespace OnlineHelp
} // namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/CoCoA-5/OnlineHelp.C,v 1.72 2014/08/01 11:23:45 bigatti Exp $
// $Log: OnlineHelp.C,v $
// Revision 1.72  2014/08/01 11:23:45  bigatti
// -- added skipping comments for ManExamples
//
// Revision 1.71  2014/07/18 12:48:07  bigatti
// -- improvement for make ManExamples
//
// Revision 1.70  2014/07/18 11:21:12  bigatti
// -- minor cleanup
//
// Revision 1.69  2014/07/18 11:08:35  bigatti
// -- removed unwanted debugging print
// -- fixed trailing spaces before comment
//
// Revision 1.68  2014/07/18 08:27:54  bigatti
// -- now "all further matches" does not print what was in seealso
//
// Revision 1.67  2014/07/16 13:24:57  bigatti
// -- added obsolete_functions and obsolescent_functions
//
// Revision 1.66  2014/05/14 07:29:37  bigatti
// -- improved (fixed) PrintExampleWithoutOutput
//
// Revision 1.65  2014/05/09 10:57:33  bigatti
// -- added -- before printing title so that it gets coloured as comments
//
// Revision 1.64  2014/04/24 11:37:34  abbott
// Summary: Added removal of <rtn> and </rtn> to CleanupLine
// Author: JAA
//
// Revision 1.63  2014/04/24 11:03:28  bigatti
// -- rtntype cahnged into rtn
//
// Revision 1.62  2014/04/23 16:48:44  bigatti
// -- fixed RtnTypes
//
// Revision 1.61  2014/04/23 14:15:38  bigatti
// ++ added myAddRetType
//
// Revision 1.60  2014/04/22 17:18:21  bigatti
// ++ removed backslash (no longer needed)
//
// Revision 1.59  2014/04/22 12:46:07  abbott
// Summary: Text inside <em>...</em> now printed between *...*
// Author: JAA
//
// Revision 1.58  2014/04/09 13:28:38  bigatti
// -- removed printing title for syntax and description
//
// Revision 1.57  2014/04/07 11:10:34  bigatti
// -- fixed matches for case ExactMatchKey==1
//
// Revision 1.56  2014/04/07 10:25:29  abbott
// Summary: Added {} as suggested by compiler; added condition to avoid duplicate in "other matches"
// Author: JAA
//
// Revision 1.55  2014/04/04 13:42:00  bigatti
// -- another fix
//
// Revision 1.54  2014/04/04 12:06:24  bigatti
// -- another fix
//
// Revision 1.53  2014/04/04 11:13:28  bigatti
// -- another fix
//
// Revision 1.52  2014/04/04 11:07:44  bigatti
// -- fix
//
// Revision 1.51  2014/04/04 11:01:04  bigatti
// -- now also searching if there is a single exact match in keys
//
// Revision 1.50  2014/03/28 15:23:41  bigatti
// -- removed "syntax" for "commands and functions for .."
//
// Revision 1.49  2014/03/27 17:33:58  bigatti
// -- added "No entry for "keyword"" message
//
// Revision 1.48  2014/03/26 16:28:04  abbott
// Summary: Improved messages printed out by online help
// Author: JAA
//
// Revision 1.47  2014/03/26 15:20:26  bigatti
// -- not printing "all matches" automatically if there are too many
//
// Revision 1.46  2014/03/26 12:02:01  abbott
// Summary: Major change to PrintMan: if search string starts with '?' then just list titles of matching manual entries
// Author: JAA
//
// Revision 1.45  2014/03/19 15:51:54  abbott
// Summary: Added blank line after syntax section; changed "See Also" --> "see also"
// Author: JAA
//
// Revision 1.44  2014/03/06 15:52:52  abbott
// Summary: Cleaned impl of PrintMan
// Author: JAA
//
// Revision 1.43  2014/03/04 14:27:54  bigatti
// -- added CleanupLine for myAddKeys
//
// Revision 1.42  2013/06/12 08:55:21  bigatti
// -- fixed "syntax" for error in entry "CartesianProduct"
//
// Revision 1.41  2013/02/27 10:44:24  bigatti
// -- added <tt> in CleanupLine (actually useless)
//
// Revision 1.40  2013/02/26 14:18:01  bigatti
// -- added field for syntax in "index"
//
// Revision 1.39  2013/02/26 13:40:31  bigatti
// -- added CleanupLine to "title"  (for the variable "it")
//
// Revision 1.38  2013/02/22 18:09:25  bigatti
// -- removed unused arg out from SkipComment
//
// Revision 1.37  2013/02/22 11:03:41  bigatti
// -- improved SkipComment: more robust
//
// Revision 1.36  2013/02/22 10:34:28  bigatti
// ++ added SkipComment (for multiline comments ONLY!)
//
// Revision 1.35  2013/02/19 18:52:50  abbott
// Added line to ignore HTML italic indications <i>...</i>.
//
// Revision 1.34  2012/06/19 15:28:47  bigatti
// -- aesthetics: changed delimiters for <example>
//
// Revision 1.33  2012/06/18 10:02:44  bigatti
// -- added mechanism to get version number and date from the xml file
//
// Revision 1.32  2012/06/04 09:34:05  bigatti
// -- added PrintWordlist
//
// Revision 1.31  2012/04/04 13:56:35  bigatti
// -- added PrintAllExamplesWithoutOutput
//
// Revision 1.30  2012/04/02 15:15:37  bigatti
// -- added check for trailing tab
//
// Revision 1.29  2012/03/20 14:34:09  bigatti
// -- added recognition of some xml tags
// -- minor cleaning
// -- Printing "matches for .." only when meaningful
//
// Revision 1.28  2012/03/13 15:59:12  bigatti
// -- xml source moved into CoCoAManual/
//
// Revision 1.27  2012/03/09 17:01:27  bigatti
// -- removed warning about old manual
// -- added "end example"
//
// Revision 1.26  2012/03/07 18:10:40  bigatti
// -- fixed interpretation of <backslash/>
//
// Revision 1.25  2012/02/24 13:10:12  bigatti
// -- added ReloadMan
//
// Revision 1.24  2012/01/25 14:33:15  bigatti
// -- added code for checking proper sorting (inactive)
//
// Revision 1.23  2011/11/02 14:41:16  bigatti
// -- added syntax
// -- removing trailing comments and spaces to searching string
// -- minor aesthetics change
//
// Revision 1.22  2011/07/06 15:52:13  bigatti
// -- improved IsSubstring
//
// Revision 1.21  2011/05/25 12:19:15  bigatti
// -- added case <par />
//
// Revision 1.20  2011/05/24 08:17:16  bigatti
// -- added warning also at the end of description
//
// Revision 1.19  2011/05/23 13:30:37  bigatti
// -- chencged include CoCoA/library --> CoCoA/error
//
// Revision 1.18  2011/05/13 10:26:37  bigatti
// -- added warning "thi is CoCoA-4 manual"
//
// Revision 1.17  2011/05/04 12:03:50  lagorio
// *** empty log message ***
//
// Revision 1.16  2011/03/23 17:31:22  bigatti
// -- added <code>..</code>
//
// Revision 1.15  2011/02/16 17:22:40  bigatti
// -- added "See Also"
//
// Revision 1.14  2011/02/16 16:14:09  bigatti
// -- class deefinition moved into .C file
// -- added storing of types, and function for <commands_and_functions_for>
// -- cleaning up
//
// Revision 1.13  2011/02/14 10:10:08  bigatti
// -- fixed &amp; &apos; <backslash/>  ...
// -- added function PrintAllExamples
//
// Revision 1.12  2011/01/25 15:45:23  bigatti
// -- added tags verbatim and formula
//
// Revision 1.11  2010/12/26 13:14:44  abbott
// Changed GlobalLogput to clog.
//
// Revision 1.10  2010/11/22 17:45:02  bigatti
// -- fixed bug in storing keys
//
// Revision 1.9  2010/11/22 17:37:44  bigatti
// -- improved message for entries with no matching
//
// Revision 1.8  2010/09/02 13:04:08  bigatti
// -- added conversion for <em>
//
// Revision 1.7  2010/09/01 13:24:48  bigatti
// -- moved all manual functions into CoCoA::OnlineHelp namespace
//
// Revision 1.6  2010/09/01 12:27:52  lagorio
// *** empty log message ***
//
// Revision 1.5  2010/09/01 09:23:19  bigatti
// -- added some "const"
// -- using CoCoA::GlobalLogput and CoCoA_ERROR instead of clog
//
// Revision 1.4  2010/09/01 08:28:35  lagorio
// *** empty log message ***
//
// Revision 1.3  2010/09/01 07:46:19  lagorio
// *** empty log message ***
//
// Revision 1.2  2010/08/31 14:55:02  bigatti
// -- pretty basic fixes
//

