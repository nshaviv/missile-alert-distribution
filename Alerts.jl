using HTTP
using JSON
using Dates, TimeZones
using DataFrames, JLD2, JSON
using JLD2


function get_json(url)
    response = HTTP.get(url)
    json = JSON.parse(String(response.body))
    return json
end

function write_jsonobj_to_jld2!(file, obj, current_path::String = "")
    for (k, v) in pairs(obj)
        path = isempty(current_path) ? string(k) : "$current_path/$k"

        if v isa AbstractDict || typeof(v).name.name == :Object
            write_jsonobj_to_jld2!(file, v, path)
        else

            file[path] = v
        end
    end
    return
end

jldopen("alerts.jld2", "w") do file
    println("GET cities and areas...")
    regions = get_json("https://www.tzevaadom.co.il/static/cities.json?v=10")

    cities = DataFrame(values(regions["cities"]))
    sort!(cities, :id)
    file["cities"] = cities

    areas = DataFrame(values(regions["areas"]))
    areas.id .= parse.(Int, keys(regions["areas"]))
    sort!(areas, :id)
    file["areas"] = areas

    println("GET polygons...")
    polygons = get_json("https://www.tzevaadom.co.il/static/polygons.json?v=10")
    write_jsonobj_to_jld2!(file, polygons, "polygons")

    println("GET alerts and parse them...")
    url = "https://www.tzevaadom.co.il/static/historical/all.json"
    response = HTTP.get(url)
    alerts = JSON.parse(String(response.body))

    df = DataFrame(alert_id = [], threat = [], city = [], date = [])
    for row in alerts
        for city in row[3]
            new_row = Dict(
                "alert_id" => Int(row[1]),
                "threat" => Int(row[2]),
                "city" => city,
                "date" => unix2datetime(row[4])
            )
            push!(df, new_row)
        end
    end
    sort!(df, :date)
    zdt_utc = ZonedDateTime.(df.date, tz"UTC")
    df.date = astimezone.(zdt_utc, tz"Asia/Jerusalem")
    file["alerts"] = df
end