{
    gSystem->Load("libG2PSim.so");

    gSystem->AddIncludePath("-I$PWD/include");
    gInterpreter->AddIncludePath("$PWD/include");
}
