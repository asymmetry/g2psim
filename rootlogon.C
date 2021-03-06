{
    gSystem->Load("libG2PSim.so");

    gSystem->AddIncludePath("-I$PWD/include");
    gInterpreter->AddIncludePath("$PWD/include");
    gSystem->AddIncludePath("-I$PWD/HRSTrans");
    gInterpreter->AddIncludePath("$PWD/HRSTrans");
    gSystem->AddIncludePath("-I$PWD/G2PPhys");
    gInterpreter->AddIncludePath("$PWD/G2PPhys");

    if (TString(gSystem->Getenv("LIBCONFIG")) != "/usr")
    {
        gSystem->AddIncludePath("-I$LIBCONFIG/include");
        gInterpreter->AddIncludePath("$LIBCONFIG/include");
    }
}
