namespace Outer { namespace Inner { int i = 0; }}

int main(int argc, char* argv[])
{
    using namespace Outer::Inner;
    return i;
}
	
