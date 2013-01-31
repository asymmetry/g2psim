{
    gSystem->Load("libG2PSim.so");

    gSystem->AddIncludePath("-I$PWD/include");
    gInterpreter->AddIncludePath("$PWD/include");
    gSystem->AddIncludePath("-I$PWD/HRSTransport");
    gInterpreter->AddIncludePath("$PWD/HRSTransport");
    gSystem->AddIncludePath("-I$PWD/G2PXSection");
    gInterpreter->AddIncludePath("$PWD/G2PXSection");
}
