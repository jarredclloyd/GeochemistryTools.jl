using Dates

@testset "automatic date-time parsing" begin
    @test automatic_datetime("2025-02-11 21:45:54") == dateformat"Y-m-d H:M:S"
    @test automatic_datetime("2025/02/11 21:45:54") == dateformat"Y/m/d H:M:S"
    @test automatic_datetime("11/02/2025 21:45:54") == dateformat"d/m/Y H:M:S"
    @test automatic_datetime("02/11/2025 21:45:54"; day_first=false) == dateformat"m/d/Y H:M:S"
    @test automatic_datetime("2025-02-11 9:45:54 PM") == dateformat"Y-m-d H:M:S p"
    @test automatic_datetime("2025/02/11 9:45:54 PM") == dateformat"Y/m/d H:M:S p"
    @test automatic_datetime("11/02/2025 9:45:54 PM") == dateformat"d/m/Y H:M:S p"
    @test automatic_datetime("02/11/2025 9:45:54 PM"; day_first=false) == dateformat"m/d/Y H:M:S p"
end
