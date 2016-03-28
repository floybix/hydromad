##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


eventsExplorer <-
    function(Q, P = NULL,
             ..., type = list(P = "s", "l"), xlab = NULL,
             xscale.components = xscale.components.subticks)
{
    if(!require("playwith")) stop("package playwith is required for eventsExplorer")
    title <- paste("eventsExplorer:", toString(deparse(substitute(Q)), width = 33))
    DATA <- Q
    if (!is.null(P))
        DATA <- merge(Q = Q, P = P)
    xyplotObj <-
        xyplot(DATA, ..., type = type, xlab = xlab, xscale.components = xscale.components)

    layer.eventseq.QP <- 
        function(Q, P = NULL,
                 ...,
                 col = c("grey90", "grey80"), border = "grey80",
                 P.nperyear = 0, Q.nperyear = 0,
                 P.qtile = 0, Q.qtile = 0,
                 in.P.qtile = 0, in.Q.qtile = 0,
                 continue = FALSE, all = FALSE,
                 mingap = 1, mindur = 1, indur = 1, extend = 0)
        {
            P.thresh <- 0
            Q.thresh <- 0
            ## calculate thresholds as quantiles
            if (P.qtile > 0)
                P.thresh <- quantile(coredata(P), P.qtile, na.rm = TRUE, names = FALSE)
            if (Q.qtile > 0)
                Q.thresh <- quantile(coredata(Q), Q.qtile, na.rm = TRUE, names = FALSE)
            in.P.thresh <- 0
            in.Q.thresh <- 0
            inthresh <- NULL
            inx <- NULL
            if (in.P.qtile > 0) {
                in.P.thresh <- inthresh <-
                    quantile(coredata(P), in.P.qtile, na.rm = TRUE, names = FALSE)
                inx <- P
            }
            if (in.Q.qtile > 0) {
                in.Q.thresh <- inthresh <-
                    quantile(coredata(Q), in.Q.qtile, na.rm = TRUE, names = FALSE)
                inx <- Q
            }
            if ((in.P.qtile > 0) && (in.Q.qtile > 0)) {
                inthresh <- c(in.Q.thresh, in.P.thresh)
                inx <- merge(Q, P)
            }
            ## calculate thresholds from target number per year
            nyears <- as.double(end(Q) - start(Q), units = "days") / 365.25
            use.P <- ((P.thresh > 0) || (P.nperyear > 0))
            if (use.P) {
                ev <- eventseq(P, thresh = P.thresh, n = round(P.nperyear * nyears),
                               inthresh = inthresh, inx = inx,
                               mingap = mingap, mindur = mindur, indur = indur,
                               extend = extend, continue = continue, all = all)
                P.thresh <- attr(ev, "thresh")
            } else {
                ev <- eventseq(Q, thresh = Q.thresh, n = round(Q.nperyear * nyears),
                               inthresh = inthresh, inx = inx,
                               mingap = mingap, mindur = mindur, indur = indur,
                               extend = extend, continue = continue, all = all)
                Q.thresh <- attr(ev, "thresh")
            }
            note <-
                paste(c(if (use.P) c(" thresh(P) = ", signif(P.thresh, 4))
                        else c(" thresh(Q) = ", signif(Q.thresh, 4)),
                        if (in.P.qtile > 0) c("\n inthresh(P) = ", signif(in.P.thresh, 4)),
                        if (in.Q.qtile > 0) c("\n inthresh(Q) = ", signif(in.Q.thresh, 4)),
                        "\n N = ", nlevels(ev), " (", round(nlevels(ev) / nyears, 1), " / year)",
                        if (!continue && !all) c("\n coverage ", round(mean(!is.na(ev))*100), "%")),
                      collapse = "")
            ## construct composite layer
            layer_(do.call("panel.xblocks", Args),
                   data = list(Args = list(ev, col = col, border = border, ...))) +
            layer({
                if (packet.number() == 1) { ## 'Q' panel
                    if (use.P == FALSE)
                        panel.abline(h = Q.thresh, lwd = 0.5)
                    if (in.Q.qtile > 0)
                        panel.abline(h = in.Q.thresh, lwd = 0.5, lty = 3)
                }
                if (packet.number() == 2) { ## 'P' panel
                    if (use.P)
                        panel.abline(h = P.thresh, lwd = 0.5)
                    if (in.P.qtile > 0)
                        panel.abline(h = in.P.thresh, lwd = 0.5, lty = 3)
                }
            },
                  data = environment()) +
            layer(panel.usertext(current.panel.limits()$x[1],
                                 current.panel.limits()$y[2],
                                 note, cex = 0.7, adj = c(0,1.2)),
                  data = list(note = note), packets = 1)
        }
    
    playwith(update(xyplotObj) +
             layer.eventseq.QP(Q, P,
                               P.nperyear = P.nperyear,
                               Q.nperyear = Q.nperyear,
                               P.qtile = P.qtile,
                               Q.qtile = Q.qtile,
                               in.P.qtile = in.P.qtile,
                               in.Q.qtile = in.Q.qtile,
                               continue = continue,
                               all = all,
                               mingap = mingap,
                               mindur = mindur,
                               indur = indur,
                               extend = extend),
             time.mode = TRUE,
             title = title,
             parameters = list(
                 P.nperyear = 0:99,
                 Q.nperyear = list(0:99, label = " (or) Q.nperyear"),
                 P.qtile = list(seq(0, 0.999, by = 0.001), label = " (or) P qtile"),
                 Q.qtile = list(seq(0, 0.999, by = 0.001), label = " (or) Q qtile"),
                 in.P.qtile = list(seq(0, 0.999, by = 0.001), label = "inner P qtile"),
                 in.Q.qtile = list(seq(0, 0.999, by = 0.001), label = "inner Q qtile"),
                 continue = FALSE,
                 all = FALSE,
                 mingap = 1:999,
                 indur = 1:999,
                 mindur = 1:999,
                 extend = 0:999)
    )
}
