# FarmaClip — Clippy-inspired analysis coach for GS → Trait DR → BLUP(PC)

clippyCoachUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML(paste0("
      #", ns("clippy_root"), " {
        position: fixed; right: 18px; bottom: 18px; z-index: 9999;
        font-family: 'Segoe UI', Tahoma, sans-serif;
      }
      #", ns("clippy_root"), " .fc-panel {
        width: 340px; max-height: 70vh; overflow: auto;
        background: #fffef5; border: 2px solid #2c3e50; border-radius: 12px;
        box-shadow: 0 8px 24px rgba(0,0,0,.18); padding: 12px 14px 10px;
      }
      #", ns("clippy_root"), " .fc-head {
        display: flex; align-items: center; gap: 10px; margin-bottom: 8px;
      }
      #", ns("clippy_root"), " .fc-face {
        width: 52px; height: 52px; border-radius: 50%;
        background: linear-gradient(145deg, #ffeaa7, #fdcb6e);
        border: 2px solid #2c3e50; display: flex; align-items: center;
        justify-content: center; font-size: 28px; flex-shrink: 0;
        animation: fc-bob 2.4s ease-in-out infinite;
      }
      @keyframes fc-bob {
        0%,100% { transform: translateY(0); }
        50% { transform: translateY(-4px); }
      }
      #", ns("clippy_root"), " .fc-title { font-weight: 700; color: #2c3e50; margin: 0; }
      #", ns("clippy_root"), " .fc-bubble {
        background: #fff; border: 1px solid #bdc3c7; border-radius: 8px;
        padding: 8px 10px; margin: 6px 0; font-size: 13px; line-height: 1.35;
      }
      #", ns("clippy_root"), " .fc-next {
        background: #eaf6ff; border-color: #3498db; font-weight: 600;
      }
      #", ns("clippy_root"), " .fc-steps { font-size: 12px; margin: 6px 0 4px; }
      #", ns("clippy_root"), " .fc-steps span.done { color: #27ae60; font-weight: 700; }
      #", ns("clippy_root"), " .fc-steps span.todo { color: #7f8c8d; }
      #", ns("clippy_root"), " .fc-mini { font-size: 11px; color: #555; }
      #", ns("clippy_root"), " .fc-toggle {
        margin-top: 6px; width: 100%;
      }
      #", ns("clippy_root"), ".fc-collapsed .fc-panel { display: none; }
      #", ns("clippy_root"), " .fc-fab {
        width: 56px; height: 56px; border-radius: 50%; border: 2px solid #2c3e50;
        background: #fdcb6e; font-size: 28px; cursor: pointer;
        box-shadow: 0 4px 12px rgba(0,0,0,.2);
      }
    "))),
    div(
      id = ns("clippy_root"),
      actionButton(ns("toggle"), label = NULL, class = "fc-fab",
                   icon = NULL, title = "FarmaClip coach"),
      uiOutput(ns("panel"))
    ),
    tags$script(HTML(sprintf("
      $(document).on('shiny:connected', function() {
        var root = document.getElementById('%s');
        var btn = document.getElementById('%s');
        if (btn) btn.innerHTML = '📎';
      });
    ", ns("clippy_root"), ns("toggle"))))
  )
}

#' @param plant_summaries reactive (or NULL)
#' @param dr_result reactive from Trait DR (or NULL)
#' @param blup_result reactive from BLUP (or NULL)
clippyCoachServer <- function(id, plant_summaries, dr_result, blup_result) {
  moduleServer(id, function(input, output, session) {
    collapsed <- reactiveVal(FALSE)

    observeEvent(input$toggle, {
      collapsed(!isTRUE(collapsed()))
    })

    status <- reactive({
      ps <- tryCatch(plant_summaries(), error = function(e) NULL)
      dr <- tryCatch(dr_result(), error = function(e) NULL)
      bl <- tryCatch(blup_result(), error = function(e) NULL)
      .sf_coach_status(ps, dr, bl, alpha = 0.05)
    })

    observe({
      session$sendCustomMessage(
        "farmaClipCollapse",
        list(id = session$ns("clippy_root"), collapsed = isTRUE(collapsed()))
      )
    })

    # collapse via class on root
    output$panel <- renderUI({
      if (isTRUE(collapsed())) return(NULL)
      st <- status()
      sd <- st$steps_done
      step_html <- tagList(
        tags$span(class = if (isTRUE(sd[["gs"]])) "done" else "todo", "1.Growth"),
        " → ",
        tags$span(class = if (isTRUE(sd[["dr"]])) "done" else "todo", "2.Trait DR"),
        " → ",
        tags$span(class = if (isTRUE(sd[["blup"]])) "done" else "todo", "3.BLUP"),
        " → ",
        tags$span(class = if (st$step >= 4L) "done" else "todo", "4.Answer")
      )

      pcs_ui <- NULL
      if (length(st$pcs)) {
        pcs_ui <- tags$div(
          class = "fc-bubble",
          tags$strong("Dimensions"),
          tags$ul(
            style = "margin:4px 0 0 16px;padding:0;",
            if (!is.null(st$pcs$PC1)) tags$li(tags$strong("PC1: "), st$pcs$PC1),
            if (!is.null(st$pcs$PC2)) tags$li(tags$strong("PC2: "), st$pcs$PC2),
            if (!is.null(st$pcs$var_note)) tags$li(st$pcs$var_note)
          )
        )
      }

      answer_ui <- NULL
      biology_ui <- NULL
      if (st$step >= 4L) {
        sig <- st$sig
        if (is.data.frame(sig) && nrow(sig) > 0L) {
          lab <- intersect(c("contrast", "cultivar", "level"), names(sig))[1]
          pcol <- intersect(c("p.value", "p.adj"), names(sig))[1]
          rows <- apply(sig, 1, function(r) {
            paste0(
              if (!is.na(lab)) r[[lab]] else "?",
              if (!is.na(pcol)) sprintf(" (p=%.3g)", as.numeric(r[[pcol]])) else ""
            )
          })
          answer_ui <- tags$div(
            class = "fc-bubble",
            style = "border-color:#27ae60;",
            tags$strong("Significant vs check"),
            tags$ul(lapply(utils::head(rows, 12L), tags$li))
          )
        } else if (is.data.frame(st$ranks_head) && nrow(st$ranks_head) > 0L) {
          top <- utils::head(st$ranks_head, 5L)
          answer_ui <- tags$div(
            class = "fc-bubble",
            tags$strong(sprintf("Top BLUP ranks on %s", st$trait_used %||% "trait")),
            tags$ul(lapply(seq_len(nrow(top)), function(i) {
              tags$li(sprintf(
                "#%s %s (BLUP=%s)",
                top$rank[i], top$level[i], signif(top$blup[i], 3)
              ))
            })),
            tags$p(class = "fc-mini", "No FDR-significant Dunnett contrasts at α=0.05 for the current check.")
          )
        }
        if (!is.null(st$biology) && length(st$biology$bullets)) {
          biology_ui <- tags$div(
            class = "fc-bubble",
            style = "border-color:#8e44ad;",
            tags$strong("Biology"),
            tags$ul(
              style = "margin:4px 0 0 16px;padding:0;font-size:12px;",
              lapply(utils::head(st$biology$bullets, 10L), tags$li)
            )
          )
        }
      }

      div(
        class = "fc-panel",
        div(
          class = "fc-head",
          div(class = "fc-face", "📎"),
          div(
            tags$p(class = "fc-title", "FarmaClip"),
            tags$p(class = "fc-mini", "It looks like you're analyzing growth curves…")
          )
        ),
        div(class = "fc-steps", step_html),
        div(class = "fc-bubble", st$message),
        div(class = "fc-bubble fc-next", tags$strong("Next: "), st$next_action),
        pcs_ui,
        answer_ui,
        biology_ui,
        tags$p(
          class = "fc-mini",
          "Left-bar cultivar/treatment filters apply to DR/BLUP. Tab controls choose the PC trait & checks."
        )
      )
    })

    # Toggle collapsed class via JS
    session$onFlushed(function() {
      insertUI(
        selector = "body",
        where = "beforeEnd",
        ui = tags$script(HTML(sprintf("
          Shiny.addCustomMessageHandler('farmaClipCollapse', function(msg) {
            var el = document.getElementById(msg.id);
            if (!el) return;
            if (msg.collapsed) el.classList.add('fc-collapsed');
            else el.classList.remove('fc-collapsed');
          });
        "))),
        immediate = TRUE
      )
    }, once = TRUE)

    invisible(list(status = status))
  })
}
