value_temp = std::vector<std::vector<double>> ((NY+NS),std::vector<double>(singleTableSize_));

count = 0;    //counter

for (int ii = 0; ii < NZ; ++ii)
{
    for (int jj = 0; jj < NC; ++jj)
    {
        for (int kk = 0; kk < NGZ; ++kk)
        {   
            for (int ll = 0; ll < NGC; ++ll)
            {
                if (std::getline(table, line))
                {
                    std::istringstream iss(line);
                    for (int ivar=0; ivar<(NY+NS); ++ivar)
                    {
                        iss >> value_temp[ivar][count];
                    }
                }
                else
                {
                    Info << "Error: no more line" << endl;
                }
                count++;
            }
        }
    }
    prog = 100.0*ii / ((NZ-1));
    Info << " Progress -- " << prog << "%" << endl;
}