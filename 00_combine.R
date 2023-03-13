# merge the title page and manuscript

qpdf::pdf_combine(
  input = c("01_title_page.pdf","01_manuscript.pdf"),
  output = "01_full_manuscript.pdf"
)
