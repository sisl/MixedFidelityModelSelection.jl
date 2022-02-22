# r(θ) = 1 + cos(θ) * sin(θ)^2
# plot(r, 0, 2π, lims=(0,1.5), proj=:polar, label=false, xticks=([0, π], ["a","b"]))

function Plots.gr_polaraxes(rmin::Real, rmax::Real, sp::Plots.Subplot)
    Plots.GR.savestate()
    xaxis = sp[:xaxis]
    yaxis = sp[:yaxis]

    α = 0:45:315
    a = α .+ 90
    sinf = sind.(a)
    cosf = cosd.(a)
    θtick_values, θtick_labels = Plots.get_ticks(sp, xaxis, update = false)
    rtick_values, rtick_labels = Plots.get_ticks(sp, yaxis, update = false)

    #draw angular grid
    if xaxis[:grid]
        Plots.gr_set_line(
            xaxis[:gridlinewidth],
            xaxis[:gridstyle],
            xaxis[:foreground_color_grid],
            sp,
        )
        Plots.gr_set_transparency(xaxis[:foreground_color_grid], xaxis[:gridalpha])
        for i in eachindex(α)
            Plots.GR.polyline([sinf[i], 0], [cosf[i], 0])
        end
    end

    #draw radial grid
    if yaxis[:grid]
        Plots.gr_set_line(
            yaxis[:gridlinewidth],
            yaxis[:gridstyle],
            yaxis[:foreground_color_grid],
            sp,
        )
        Plots.gr_set_transparency(yaxis[:foreground_color_grid], yaxis[:gridalpha])
        for i in eachindex(rtick_values)
            r = (rtick_values[i] - rmin) / (rmax - rmin)
            if r <= 1.0 && r >= 0.0
                Plots.GR.drawarc(-r, r, -r, r, 0, 359)
            end
        end
        Plots.GR.drawarc(-1, 1, -1, 1, 0, 359)
    end

    #prepare to draw ticks
    Plots.gr_set_transparency(1)
    Plots.GR.setlinecolorind(90)
    Plots.GR.settextalign(GR.TEXT_HALIGN_CENTER, GR.TEXT_VALIGN_HALF)

    #draw angular ticks
    if xaxis[:showaxis]
        Plots.GR.drawarc(-1, 1, -1, 1, 0, 359)
        # for i in eachindex(α)
        for i in eachindex(θtick_values)
            θ = θtick_values[i]
            x, y = Plots.GR.wctondc(1.1 * sin(θ), 1.1 * cos(θ))
            # x, y = Plots.GR.wctondc(1.1 * sinf[i], 1.1 * cosf[i])
            # Plots.GR.textext(x, y, string((360 - α[i]) % 360, "^o"))
            Plots.GR.textext(x, y, Plots._cycle(θtick_labels, i))
        end
    end

    #draw radial ticks
    if yaxis[:showaxis]
        for i in eachindex(rtick_values)
            r = (rtick_values[i] - rmin) / (rmax - rmin)
            if r <= 1.0 && r >= 0.0
                x, y = Plots.GR.wctondc(0.05, r)
                Plots.gr_text(x, y, Plots._cycle(rtick_labels, i))
            end
        end
    end
    @info "Moss!"
    Plots.GR.restorestate()
end