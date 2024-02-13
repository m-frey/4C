/*----------------------------------------------------------------------*/
/*! \file
\brief Unittests for the Pstream and Level classes
\level 1

*-----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "baci_io_pstream.hpp"

#include <Epetra_SerialComm.h>

#include <stdexcept>

namespace
{
  using namespace BACI;

  TEST(PstreamTest, UninitializedUseThrows)
  {
    IO::Pstream ps;
    EXPECT_THROW(ps.flush(), CORE::Exception);
    EXPECT_THROW((ps << "blub"), CORE::Exception);
    EXPECT_NO_THROW(ps.close());
  }

  TEST(PstreamTest, DoubleInitializeThrows)
  {
    IO::Pstream ps;
    ps.setup(true, false, true, IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 4, "");
    EXPECT_THROW(
        ps.setup(false, false, false, IO::standard, Teuchos::rcp(new Epetra_SerialComm), 0, 2, ""),
        CORE::Exception);
  }

  TEST(PstreamTest, NonexistantProc)
  {
    using namespace BACI;
    IO::Pstream ps;
    EXPECT_THROW(
        ps.setup(false, false, false, IO::standard, Teuchos::rcp(new Epetra_SerialComm), 4, 2, ""),
        CORE::Exception);
  }

  TEST(PstreamTest, InitializedUse)
  {
    using namespace BACI;
    IO::Pstream ps;
    ps.setup(true, false, false, IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    EXPECT_NO_THROW(ps.flush());
    EXPECT_NO_THROW(ps << "blub");
    EXPECT_NO_THROW(ps.close());
  }

  TEST(PstreamTest, OutputLevel)
  {
    using namespace BACI;
    IO::Pstream ps;
    ps.setup(true, false, false, IO::minimal, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    EXPECT_EQ(ps.RequestedOutputLevel(), IO::minimal);
    IO::Level &lvl = ps(IO::debug);
    EXPECT_NO_THROW(lvl << 4);
  }

  TEST(PstreamTest, InputTypes)
  {
    using namespace BACI;
    IO::Pstream ps;
    ps.setup(false, false, true, IO::debug, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    EXPECT_NO_THROW(ps << 4UL << -5LL << 1337.0 << 42.0f << "blub" << std::string("blah") << "\n");
    EXPECT_NO_THROW(ps.flush());
    EXPECT_NO_THROW(ps.close());
  }

  TEST(PstreamTest, ExternalOperators)
  {
    using namespace BACI;
    IO::Pstream ps;
    ps.setup(false, false, true, IO::debug, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    EXPECT_NO_THROW(ps << "blub" << IO::flush);
    EXPECT_NO_THROW(ps << "blah" << IO::endl);
  }

  TEST(PstreamTest, Level)
  {
    using namespace BACI;
    IO::Pstream ps;
    ps.setup(true, false, true, IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    IO::Level &lvl = ps(IO::debug);
    EXPECT_NO_THROW(lvl.stream(1.2));
    EXPECT_NO_THROW(lvl << 4);
    EXPECT_NO_THROW(lvl.SetLevel(IO::minimal) << 5);
  }

  TEST(PstreamTest, LevelExternalOperators)
  {
    using namespace BACI;
    IO::Pstream ps;
    ps.setup(true, false, true, IO::standard, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    IO::Level &lvl = ps(IO::debug);
    EXPECT_NO_THROW(lvl << 1.2 << IO::flush);
    EXPECT_NO_THROW(lvl << 23 << IO::endl);
  }

}  // namespace