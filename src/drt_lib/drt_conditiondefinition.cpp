/*!----------------------------------------------------------------------
\file drt_conditiondefinition.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <algorithm>
#include <iterator>

#include "drt_conditiondefinition.H"
#include "drt_colors.H"
#include "drt_globalproblem.H"
#include "drt_discret.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::ConditionComponent::ConditionComponent(std::string name)
  : name_(name)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::ConditionComponent::PushBack(std::string token, Teuchos::RCP<std::stringstream> stream)
{
  Teuchos::RCP<std::stringstream> out = Teuchos::rcp(new std::stringstream());
  (*out) << token << " ";
  std::copy(std::istream_iterator<std::string>(*stream),
            std::istream_iterator<std::string>(),
            std::ostream_iterator<std::string>(*out," "));
  return out;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::StringConditionComponent::StringConditionComponent(std::string name,
                                                               std::string defaultvalue,
                                                               const Teuchos::Array<std::string>& datfilevalues,
                                                               const Teuchos::Array<std::string>& condvalues,
                                                               bool optional)
  : ConditionComponent(name),
    defaultvalue_(defaultvalue),
    datfilevalues_(datfilevalues),
    condvalues_(condvalues),
    optional_(optional)
{
  if (std::find(datfilevalues_.begin(),datfilevalues_.end(),defaultvalue_)==datfilevalues_.end())
  {
    dserror("invalid default value '%s'", defaultvalue_.c_str());
  }
  if (datfilevalues_.size()!=condvalues_.size())
  {
    dserror("dat file values must match condition values");
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::StringConditionComponent::DefaultLine(std::ostream& stream)
{
  stream << defaultvalue_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::StringConditionComponent::Print(std::ostream& stream, const Condition* cond)
{
  stream << *cond->Get<std::string>(Name());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::StringConditionComponent::Read(DRT::INPUT::ConditionDefinition* def,
                                                                           Teuchos::RCP<std::stringstream> condline,
                                                                           Teuchos::RCP<DRT::Condition> condition)
{
  std::string value;
  (*condline) >> value;

  if (value=="")
    value = defaultvalue_;

  Teuchos::Array<std::string>::iterator i = std::find(datfilevalues_.begin(),datfilevalues_.end(),value);
  if (i==datfilevalues_.end())
  {
    if (optional_)
    {
      condline = PushBack(value,condline);
      i = std::find(datfilevalues_.begin(),datfilevalues_.end(),defaultvalue_);
    }
    else
      dserror("unrecognized string '%s' while reading variable '%s' in '%s'",
              value.c_str(),Name().c_str(),def->SectionName().c_str());
  }
  unsigned pos = &*i - &datfilevalues_[0];
  condition->Add(Name(),condvalues_[pos]);

  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::SeparatorConditionComponent::SeparatorConditionComponent(std::string separator)
  : ConditionComponent("*SEPARATOR*"), separator_(separator)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::SeparatorConditionComponent::DefaultLine(std::ostream& stream)
{
  stream << separator_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::SeparatorConditionComponent::Print(std::ostream& stream, const DRT::Condition* cond)
{
  stream << separator_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::SeparatorConditionComponent::Read(DRT::INPUT::ConditionDefinition* def,
                                                                              Teuchos::RCP<std::stringstream> condline,
                                                                              Teuchos::RCP<DRT::Condition> condition)
{
  std::string sep;
  (*condline) >> sep;
  if (sep!=separator_)
    dserror("word '%s' expected but found '%s' while reading '%s'",
            separator_.c_str(),sep.c_str(),def->SectionName().c_str());
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::IntConditionComponent::IntConditionComponent(std::string name, bool fortranstyle, bool noneallowed)
  : ConditionComponent(name),
    fortranstyle_(fortranstyle),
    noneallowed_(noneallowed)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntConditionComponent::DefaultLine(std::ostream& stream)
{
  if (noneallowed_)
    stream << "none";
  else
    stream << 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntConditionComponent::Print(std::ostream& stream, const DRT::Condition* cond)
{
  int n = cond->GetInt(Name());
  if (noneallowed_ and n==-1)
    stream << "none ";
  else
  {
    if (fortranstyle_)
      n += 1;
    stream << n;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::IntConditionComponent::Read(DRT::INPUT::ConditionDefinition* def,
                                                                        Teuchos::RCP<std::stringstream> condline,
                                                                        Teuchos::RCP<DRT::Condition> condition)
{
  std::string number;
  (*condline) >> number;

  int n;
  if (noneallowed_ and number=="none")
  {
    n = -1;
  }
  else
  {
    char* ptr;
    n = strtol(number.c_str(),&ptr,10);
    if (ptr==number.c_str())
      dserror("failed to read number '%s' while reading variable '%s' in '%s'",
              number.c_str(),Name().c_str(),def->SectionName().c_str());
  }
  if (fortranstyle_)
  {
    if (not noneallowed_ or n!=-1)
      n -= 1;
  }
  condition->Add(Name(),n);
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::IntVectorConditionComponent::IntVectorConditionComponent(std::string name,
                                                                     int length,
                                                                     bool fortranstyle,
                                                                     bool noneallowed,
                                                                     bool optional)
  : ConditionComponent(name),
    length_(length),
    fortranstyle_(fortranstyle),
    noneallowed_(noneallowed),
    optional_(optional)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorConditionComponent::DefaultLine(std::ostream& stream)
{
  if (noneallowed_)
  {
    for (int i=0; i<length_; ++i)
      stream << "none ";
  }
  else
  {
    for (int i=0; i<length_; ++i)
      stream << 0 << " ";
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorConditionComponent::Print(std::ostream& stream, const DRT::Condition* cond)
{
  const std::vector<int>* v = cond->Get<std::vector<int> >(Name());
  for (unsigned i=0; i<v->size(); ++i)
  {
    if (noneallowed_ and (*v)[i]==-1)
      stream << "none ";
    else
    {
      if (fortranstyle_)
        stream << (*v)[i]+1 << " ";
      else
        stream << (*v)[i]+1 << " ";
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::IntVectorConditionComponent::Read(DRT::INPUT::ConditionDefinition* def,
                                                                              Teuchos::RCP<std::stringstream> condline,
                                                                              Teuchos::RCP<DRT::Condition> condition)
{
  std::vector<int> numbers(length_,0);

  for (int i=0; i<length_; ++i)
  {
    std::string number;
    (*condline) >> number;

    int n;
    if (noneallowed_ and number=="none")
    {
      n = -1;
    }
    else
    {
      char* ptr;
      n = strtol(number.c_str(),&ptr,10);
      if (ptr==number.c_str())
      {
        if (optional_ and i==0)
        {
          // failed to read the numbers, fall back to default values
          condline = PushBack(number,condline);
          break;
        }
        dserror("failed to read number '%s' while reading variable '%s' in '%s'",
                number.c_str(),Name().c_str(),def->SectionName().c_str());
      }
    }
    if (fortranstyle_)
    {
      if (not noneallowed_ or n!=-1)
        n -= 1;
    }
    numbers[i] = n;
  }
  condition->Add(Name(),numbers);
  return condline;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorConditionComponent::SetLength(int length)
{
  length_=length;
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::RealConditionComponent::RealConditionComponent(std::string name)
  : ConditionComponent(name)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealConditionComponent::DefaultLine(std::ostream& stream)
{
  stream << "0.0";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealConditionComponent::Print(std::ostream& stream, const DRT::Condition* cond)
{
  stream << cond->GetDouble(Name());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::RealConditionComponent::Read(DRT::INPUT::ConditionDefinition* def,
                                                                         Teuchos::RCP<std::stringstream> condline,
                                                                         Teuchos::RCP<DRT::Condition> condition)
{
  double number = 0;
  (*condline) >> number;
  condition->Add(Name(),number);
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::RealVectorConditionComponent::RealVectorConditionComponent(std::string name, int length)
  : ConditionComponent(name), length_(length)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorConditionComponent::DefaultLine(std::ostream& stream)
{
  for (int i=0; i<length_; ++i)
    stream << "0.0 ";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorConditionComponent::Print(std::ostream& stream, const DRT::Condition* cond)
{
  const std::vector<double>* v = cond->Get<std::vector<double> >(Name());
  for (unsigned i=0; i<v->size(); ++i)
    stream << (*v)[i] << " ";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::RealVectorConditionComponent::Read(DRT::INPUT::ConditionDefinition* def,
                                                                               Teuchos::RCP<std::stringstream> condline,
                                                                               Teuchos::RCP<DRT::Condition> condition)
{
  std::vector<double> numbers(length_);

  for (int i=0; i<length_; ++i)
  {
    double number = 0;
    (*condline) >> number;
    numbers[i] = number;
  }
  condition->Add(Name(),numbers);
  return condline;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorConditionComponent::SetLength(int length)
{
  length_=length;
  return;
}

DRT::INPUT::DirichletNeumannBundle::DirichletNeumannBundle
(
  std::string name, 
  Teuchos::RCP<IntConditionComponent> intcomp,
  std::vector<Teuchos::RCP<IntVectorConditionComponent> > intvectcomp,
  std::vector<Teuchos::RCP<RealVectorConditionComponent> > realvectcomp
)
  : ConditionComponent(name),
    intcomp_(intcomp),
    intvectcomp_(intvectcomp),
    realvectcomp_(realvectcomp)
{};

void DRT::INPUT::DirichletNeumannBundle::DefaultLine(std::ostream& stream)
{
  intcomp_->DefaultLine(stream);
  stream << "  ";
  intvectcomp_[0]->DefaultLine(stream);
  stream << " ";
  realvectcomp_[0]->DefaultLine(stream);
  stream << " ";
  intvectcomp_[1]->DefaultLine(stream);
  stream << " ";
  intvectcomp_[2]->DefaultLine(stream);
  stream << " ";
}

void DRT::INPUT::DirichletNeumannBundle::Print(std::ostream& stream, const DRT::Condition* cond)
{
  intcomp_->Print(stream,cond);
  stream << "  ";
  intvectcomp_[0]->Print(stream,cond);
  stream << " ";
  realvectcomp_[0]->Print(stream,cond);
  stream << " ";
  intvectcomp_[1]->Print(stream,cond);
  stream << " ";
  intvectcomp_[2]->Print(stream,cond);
  stream << " ";
};

Teuchos::RCP<std::stringstream> DRT::INPUT::DirichletNeumannBundle::Read(ConditionDefinition* def,
                                             Teuchos::RCP<std::stringstream> condline,
                                             Teuchos::RCP<DRT::Condition> condition)
{ 
  intcomp_->Read(def,condline,condition);
  int length = condition->GetInt(intcomp_->Name());
  
  intvectcomp_[0]->SetLength(length);
  intvectcomp_[0]->Read(def,condline,condition);
  realvectcomp_[0]->SetLength(length);
  realvectcomp_[0]->Read(def,condline,condition);
  intvectcomp_[1]->SetLength(length);
  intvectcomp_[1]->Read(def,condline,condition);
  intvectcomp_[2]->SetLength(length);
  intvectcomp_[2]->Read(def,condline,condition);
  
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::ConditionDefinition::ConditionDefinition(std::string sectionname,
                                                     std::string conditionname,
                                                     std::string description,
                                                     DRT::Condition::ConditionType condtype,
                                                     bool buildgeometry,
                                                     DRT::Condition::GeometryType gtype)
  : sectionname_(sectionname),
    conditionname_(conditionname),
    description_(description),
    condtype_(condtype),
    buildgeometry_(buildgeometry),
    gtype_(gtype)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::ConditionDefinition::AddComponent(Teuchos::RCP<ConditionComponent> c)
{
  inputline_.push_back(c);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::ConditionDefinition::Read(const Problem& problem,
    DatFileReader& reader,
    std::multimap<int,Teuchos::RCP<DRT::Condition> >& cmap)
{
  std::string name = "--";
  name += sectionname_;
  std::vector<const char*> section = reader.Section(name);

  if (section.size()>0)
  {
    std::stringstream line(section[0]);
    std::string dobj;
    int condcount=-1;
    line >> dobj;
    line >> condcount;

    if(condcount<0)
    {
      dserror("condcount<0 in section %s\n",sectionname_.c_str());
    }

    bool success = false;
    switch (gtype_)
    {
    case DRT::Condition::Point:   success = dobj=="DPOINT"; break;
    case DRT::Condition::Line:    success = dobj=="DLINE";  break;
    case DRT::Condition::Surface: success = dobj=="DSURF";  break;
    case DRT::Condition::Volume:  success = dobj=="DVOL";   break;
    default:
      dserror("geometry type unspecified");
    }

    if (not success)
    {
      dserror("expected design object type but got '%s' in '%s'",
              dobj.c_str(),sectionname_.c_str());
    }

    if (condcount != static_cast<int>(section.size()-1))
    {
      dserror("Got %d condition lines but expect %d in '%s'",
              section.size()-1,condcount,sectionname_.c_str());
    }

    const Teuchos::ParameterList& conditionnames = problem.ConditionNamesParams();

    for (std::vector<const char*>::iterator i=section.begin()+1;
         i!=section.end();
         ++i)
    {
      Teuchos::RCP<std::stringstream> condline = Teuchos::rcp(new std::stringstream(*i));

      std::string E;
      std::string number;
      std::string minus;

      (*condline) >> E >> number >> minus;
      if (not (*condline) or E!="E" or minus!="-")
        dserror("invalid condition line in '%s'",sectionname_.c_str());

      int dobjid = -1;
      if (conditionnames.isParameter(number))
        dobjid = conditionnames.get<int>(number) - 1;
      else
      {
        char* ptr;
        dobjid = strtol(number.c_str(),&ptr,10);
        if (ptr==number.c_str())
          dserror("failed to read design object number '%s' in '%s'",
                  number.c_str(),sectionname_.c_str());
        dobjid -= 1;
      }

      //cout << "condition " << dobjid << " of type " << condtype_ << " on " << gtype_ << "\n";

      Teuchos::RCP<DRT::Condition> condition =
        Teuchos::rcp(new DRT::Condition(dobjid,condtype_,buildgeometry_,gtype_));

      for (unsigned j=0; j<inputline_.size(); ++j)
      {
        condline = inputline_[j]->Read(this,condline,condition);
      }

      //------------------------------- put condition in map of conditions
      cmap.insert(pair<int,Teuchos::RCP<DRT::Condition> >(dobjid,condition));
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& DRT::INPUT::ConditionDefinition::Print(std::ostream& stream,
                                                     const Discretization* dis,
                                                     bool color)
{
  std::string blue2light = "";
  std::string bluelight = "";
  std::string redlight = "";
  std::string yellowlight = "";
  std::string greenlight = "";
  std::string magentalight = "";
  std::string endcolor = "";

  if (color)
  {
    blue2light = BLUE2_LIGHT;
    bluelight = BLUE_LIGHT;
    redlight = RED_LIGHT;
    yellowlight = YELLOW_LIGHT;
    greenlight = GREEN_LIGHT;
    magentalight = MAGENTA_LIGHT;
    endcolor = END_COLOR;
  }

  unsigned l = sectionname_.length();
  stream << redlight << "--";
  for (int i=0; i<std::max<int>(65-l,0); ++i) stream << '-';
  stream << greenlight << sectionname_ << endcolor << '\n';

  std::string name;
  switch (gtype_)
  {
  case DRT::Condition::Point:   name = "DPOINT"; break;
  case DRT::Condition::Line:    name = "DLINE";  break;
  case DRT::Condition::Surface: name = "DSURF";  break;
  case DRT::Condition::Volume:  name = "DVOL";   break;
  default:
    dserror("geometry type unspecified");
  }

  int count = 0;
  if (dis!=NULL)
  {
    std::vector<Condition*> conds;
    dis->GetCondition(conditionname_,conds);
    for (unsigned c=0; c<conds.size(); ++c)
    {
      if (conds[c]->GType()==gtype_)
      {
        count += 1;
      }
    }
  }

  stream << bluelight << name << endcolor;
  l = name.length();
  for (int i=0; i<std::max<int>(31-l,0); ++i) stream << ' ';
  stream << ' ' << yellowlight << count << endcolor << '\n';

  stream << blue2light << "//" << magentalight << "E num - ";
  for (unsigned i=0; i<inputline_.size(); ++i)
  {
    inputline_[i]->DefaultLine(stream);
    stream << " ";
  }

  stream << endcolor << "\n";

  if (dis!=NULL)
  {
    std::vector<Condition*> conds;
    dis->GetCondition(conditionname_,conds);

    for (unsigned c=0; c<conds.size(); ++c)
    {
      if (conds[c]->GType()==gtype_)
      {
        stream << "E " << conds[c]->Id() << " - ";
        for (unsigned i=0; i<inputline_.size(); ++i)
        {
          inputline_[i]->Print(stream,conds[c]);
          stream << " ";
        }
        stream << "\n";
      }
    }
  }

  return stream;
}



#endif
